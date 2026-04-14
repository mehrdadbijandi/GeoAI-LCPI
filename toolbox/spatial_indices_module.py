import arcpy
import math
import os
import collections

# ==========================================
# 1. Geometric Helper Functions
# ==========================================

def get_angle_degrees(p1, p2, p3):
    """Calculates angle between 3 points for vertex cleaning."""
    if not (p1 and p2 and p3): return 0
    theta1 = math.atan2(p1.Y - p2.Y, p1.X - p2.X)
    theta2 = math.atan2(p3.Y - p2.Y, p3.X - p2.X)
    diff = theta2 - theta1
    while diff < 0: diff += 2 * math.pi
    return math.degrees(diff)

def clean_vertices_circular(part_obj, tolerance=2.0):
    """Removes collinear vertices to identify true corners for BC Index."""
    pts = [p for p in part_obj if p]
    if not pts: return []
    if len(pts) > 1:
        if pts[0].X != pts[-1].X or pts[0].Y != pts[-1].Y:
            pts.append(pts[0])
    changed = True
    while changed:
        if len(pts) < 4: break 
        changed = False
        new_pts = []
        n = len(pts) - 1
        indices_to_remove = set()
        for i in range(n):
            p_prev, p_curr, p_next = pts[i-1], pts[i], pts[(i+1) % n]
            try:
                angle = get_angle_degrees(p_prev, p_curr, p_next)
                if abs(angle - 180.0) < tolerance:
                    indices_to_remove.add(i)
                    changed = True
            except: pass
        if changed:
            for i in range(n):
                if i not in indices_to_remove: new_pts.append(pts[i])
            if len(new_pts) > 0: new_pts.append(new_pts[0])
            pts = new_pts
    return pts[:-1] if len(pts) > 1 else pts

def calculate_bc_index(parcel_perimeter, ideal_perimeter, n_vertices):
    """
    Calculates Boundary Complexity (BC).
    Formula: BC = (P_ideal / P_parcel) * (4 / N_vertices) 
    """
    if parcel_perimeter <= 0 or n_vertices == 0: return 0
    return (ideal_perimeter / parcel_perimeter) * (4.0 / n_vertices)

def get_ideal_dims(area, ratio_val):
    """Calculates ideal rectangle dimensions based on area and ratio."""
    if area <= 0: return 0, 0, 0
    try:
        r_parts = ratio_val.split(':')
        numeric_ratio = float(r_parts[0]) / float(r_parts[1])
    except: numeric_ratio = 0.5
    ideal_length = math.sqrt(area / numeric_ratio)
    ideal_width = ideal_length * numeric_ratio
    return ideal_length, ideal_width, 2 * (ideal_length + ideal_width)

def create_ideal_rect_geom(parcel_geom, length, width):
    """Creates a rotated ideal rectangle at parcel centroid."""
    try:
        mbr = parcel_geom.hullRectangle.split(' ')
        pts = [(float(mbr[i]), float(mbr[i+1])) for i in range(0, len(mbr), 2)]
        if len(pts) < 3: return None
        d12 = math.hypot(pts[1][0]-pts[0][0], pts[1][1]-pts[0][1])
        d23 = math.hypot(pts[2][0]-pts[1][0], pts[2][1]-pts[1][1])
        dx, dy = (pts[1][0]-pts[0][0], pts[1][1]-pts[0][1]) if d12 > d23 else (pts[2][0]-pts[1][0], pts[2][1]-pts[1][1])
        angle = math.atan2(dy, dx)
        cx, cy = parcel_geom.centroid.X, parcel_geom.centroid.Y
        cos_a, sin_a = math.cos(angle), math.sin(angle)
        hw, hl = width/2.0, length/2.0
        final_pts = []
        for x, y in [(-hl,-hw), (hl,-hw), (hl,hw), (-hl,hw), (-hl,-hw)]:
            final_pts.append(arcpy.Point(x*cos_a - y*sin_a + cx, x*sin_a + y*cos_a + cy))
        return arcpy.Polygon(arcpy.Array(final_pts), parcel_geom.spatialReference)
    except: return None

# ==========================================
# 2. AC Pre-processing
# ==========================================

def preprocess_accessibility_data(out_fc, in_roads, road_width_field, search_dist):
    """
    Pre-processes data for Accessibility Index (AC).
    1. Calculates 'B_max' (Maximum Road Width) .
    2. Buffers roads by (Width/2 + SearchDist) .
    3. Intersects parcel boundaries with buffers to get 'li' (effective frontage) and 'Bi' (road width).
    Returns:
        ac_data: Dictionary {OID: Sum(li * Bi)}
        w_max: Maximum road width
    """
    arcpy.AddMessage("    - Calculating road-parcel intersections...")
    
    # 1. Get Max Road Width (B_max)
    w_max = 10.0
    try:
        stats = arcpy.analysis.Statistics(in_roads, "in_memory\\Stats", [[road_width_field, "MAX"]])
        with arcpy.da.SearchCursor(stats, [f"MAX_{road_width_field}"]) as sc:
            for r in sc: w_max = float(r[0]) if r[0] else 10.0
    except: pass
    
    # Get OID field name
    desc = arcpy.Describe(out_fc)
    oid_field_name = desc.OIDFieldName
    
    # Temporary Link ID to track parcels after conversion to line
    link_fld = "Link_ID_Temp"
    if link_fld not in [f.name for f in arcpy.ListFields(out_fc)]:
        arcpy.management.AddField(out_fc, link_fld, "LONG")
    arcpy.management.CalculateField(out_fc, link_fld, f"!{oid_field_name}!", "PYTHON3")
    
    # 2. Create Road Buffers
    # Buffer Distance = (Width / 2) + Search Distance 
    temp_roads = arcpy.management.CopyFeatures(in_roads, "in_memory\\TRoads")
    arcpy.management.AddField(temp_roads, "BDist", "DOUBLE")
    arcpy.management.CalculateField(temp_roads, "BDist", f"(!{road_width_field}! / 2.0) + {search_dist}", "PYTHON3")
    road_bufs = arcpy.analysis.Buffer(temp_roads, "in_memory\\RBufs", "BDist")
    
    # 3. Convert Parcels to Line (to measure frontage length)
    p_lines = arcpy.management.PolygonToLine(out_fc, "in_memory\\PLines", "IGNORE_NEIGHBORS")
    
    # 4. Intersect
    inter = arcpy.analysis.Intersect([p_lines, road_bufs], "in_memory\\Inter")
    
    # Store Sum(li * Bi) for each parcel
    ac_data = collections.defaultdict(float)
    
    with arcpy.da.SearchCursor(inter, [link_fld, "SHAPE@LENGTH", road_width_field]) as cur:
        for row in cur:
            if row[0] is not None:
                li = row[1]  # Length of segment inside buffer
                Bi = float(row[2]) if row[2] else 0  # Width of the road
                ac_data[row[0]] += (li * Bi)
                
    # Cleanup
    try: arcpy.management.DeleteField(out_fc, link_fld)
    except: pass
    arcpy.management.Delete("in_memory")
    
    return ac_data, w_max

# ==========================================
# 3. Main Logic
# ==========================================

def main():
    arcpy.env.overwriteOutput = True
    
    # Inputs (Strictly matching the GUI)
    in_parcels = arcpy.GetParameterAsText(0)
    out_fc = arcpy.GetParameterAsText(1)
    owner_fld = arcpy.GetParameterAsText(2)
    in_roads = arcpy.GetParameterAsText(3)
    road_w_fld = arcpy.GetParameterAsText(4)
    ratio_str = arcpy.GetParameterAsText(5) or "1:2"
    out_rects = arcpy.GetParameterAsText(6)
    tol = float(arcpy.GetParameterAsText(7)) if arcpy.GetParameterAsText(7) else 2.0
    k_fact = float(arcpy.GetParameterAsText(8)) if arcpy.GetParameterAsText(8) else 500.0
    search_d = float(arcpy.GetParameterAsText(9)) if arcpy.GetParameterAsText(9) else 5.0

    # Initialize Output
    arcpy.management.CopyFeatures(in_parcels, out_fc)
    
    idx_flds = ["SI_Index", "BC_Index", "TF_Index", "AC_Index", "Dist_To_Ref"]
    for f in idx_flds:
        if f not in [fld.name for fld in arcpy.ListFields(out_fc)]:
            arcpy.management.AddField(out_fc, f, "DOUBLE")

    # --- AC Preparation ---
    ac_map = {}
    w_max = 10.0
    near_dist_field = "NEAR_DIST"
    has_roads = False

    if in_roads and road_w_fld:
        has_roads = True
        # 1. Get Direct Access data (Sum li*Bi)
        ac_map, w_max = preprocess_accessibility_data(out_fc, in_roads, road_w_fld, search_d)
        
        # 2. Calculate Distance to Nearest Road for 'No Access' parcels
        # This adds 'NEAR_DIST' field to out_fc
        arcpy.AddMessage("    - Calculating distances to nearest road...")
        arcpy.analysis.Near(out_fc, in_roads)

    # Data structures for TF calculation
    owner_groups = collections.defaultdict(list)
    rect_list = []
    global_max_area = 0.0
    
    # List to store AC values for parcels with direct road access
    direct_access_ac_values = []

    # Prepare Fields for Cursor
    fields = ["SHAPE@", "SI_Index", "BC_Index", "AC_Index", "OID@", owner_fld]
    if has_roads:
        fields.append(near_dist_field) # Add NEAR_DIST to read distances
    
    # First pass: Calculate AC for parcels with direct road access and store values
    with arcpy.da.UpdateCursor(out_fc, fields) as cursor:
        for row in cursor:
            geom, oid, owner = row[0], row[4], row[5]
            if not geom or geom.area <= 0: continue
            
            if geom.area > global_max_area: global_max_area = geom.area

            # --- SI Calculation ---
            i_len, i_wid, i_p = get_ideal_dims(geom.area, ratio_str)
            ideal_rect = create_ideal_rect_geom(geom, i_len, i_wid)
            si = 0
            if ideal_rect:
                if out_rects: rect_list.append((ideal_rect, oid))
                si = min(geom.intersect(ideal_rect, 4).area / geom.area, 1.0)
            
            # --- BC Calculation ---
            cleaned = clean_vertices_circular(geom.getPart(0), tol)
            bc = calculate_bc_index(geom.length, i_p, len(cleaned))
            
            # --- AC Calculation (First Pass) ---
            ac = 0.0
            if has_roads:
                # Check for Direct Access
                if oid in ac_map:
                    # Formula Eq 5: AC = (2 / (P * Bmax)) * Sum(li * Bi)
                    # ac_map[oid] holds Sum(li * Bi)
                    numerator = 2.0 * ac_map[oid]
                    denominator = geom.length * w_max
                    if denominator > 0:
                        ac = numerator / denominator
                        # Store AC value for direct access parcels
                        direct_access_ac_values.append(ac)
                else:
                    # No Direct Access - set AC to None for now
                    # Will be calculated in second pass with dynamic epsilon
                    ac = None
            
            # --- Gather Data for TF ---
            if owner is not None:
                owner_groups[owner].append({'oid': oid, 'shape': geom, 'area': geom.area})
            
            # Update row
            row[1], row[2], row[3] = si, bc, ac
            cursor.updateRow(row)

    # Calculate epsilon (minimum AC value from direct access parcels)
    if has_roads:
        if direct_access_ac_values:
            # Epsilon is the minimum AC value from parcels with direct road access
            epsilon = min(direct_access_ac_values)
            arcpy.AddMessage(f"    - Dynamic epsilon calculated: {epsilon:.6f} (from {len(direct_access_ac_values)} parcels with direct road access)")
        else:
            # No parcels with direct road access - use default value
            epsilon = 0.1
            arcpy.AddMessage(f"    - No parcels with direct road access. Using default epsilon: {epsilon}")
        
        # Second pass: Calculate AC for parcels without direct road access using dynamic epsilon
        arcpy.AddMessage("    - Calculating AC for parcels without direct road access...")
        with arcpy.da.UpdateCursor(out_fc, ["OID@", "AC_Index", near_dist_field]) as cursor:
            for row in cursor:
                oid, ac_val, d_road = row[0], row[1], row[2]
                
                # Only update parcels where AC is None (no direct road access)
                if ac_val is None:
                    if d_road is None:
                        d_road = 0
                    # AC = epsilon * (K / (K + D_road))
                    ac_val = epsilon * (k_fact / (k_fact + d_road))
                    row[1] = ac_val
                    cursor.updateRow(row)

    # Save Ideal Rectangles
    if out_rects and rect_list:
        path, name = os.path.split(out_rects)
        if arcpy.Exists(out_rects): arcpy.management.Delete(out_rects)
        arcpy.management.CreateFeatureclass(path, name, "POLYGON", spatial_reference=arcpy.Describe(in_parcels).spatialReference)
        arcpy.management.AddField(out_rects, "Ref_Parcel_OID", "LONG")
        with arcpy.da.InsertCursor(out_rects, ["SHAPE@", "Ref_Parcel_OID"]) as ic:
            for r in rect_list: ic.insertRow(r)

    # ==========================================
    # 4. TF Calculation 
    # ==========================================
    # TF = (A_R / A_max) * Sum( K / (K + d_j) ) 
    
    tf_results = {}
    dist_results = {}

    for owner, parcels in owner_groups.items():
        if not parcels: continue

        # Reference Parcel (R): Largest parcel of the owner
        ref_parcel = max(parcels, key=lambda x: x['area'])
        A_R = ref_parcel['area']
        ref_centroid = ref_parcel['shape'].centroid
        
        sum_distance_score = 0.0
        
        for p in parcels:
            d = 0.0
            # If parcel is reference or adjacent/touching, d = 0
            if p['oid'] == ref_parcel['oid']:
                d = 0.0
            elif not p['shape'].disjoint(ref_parcel['shape']):
                d = 0.0
            else:
                p_c = p['shape'].centroid
                d = math.hypot(p_c.X - ref_centroid.X, p_c.Y - ref_centroid.Y)
            
            dist_results[p['oid']] = d
            sum_distance_score += (k_fact / (k_fact + d))

        # TF Calculation
        final_tf = 0
        if global_max_area > 0:
            final_tf = (A_R / global_max_area) * sum_distance_score

        for p in parcels:
            tf_results[p['oid']] = final_tf

    # Update TF and Distance fields
    with arcpy.da.UpdateCursor(out_fc, ["OID@", "TF_Index", "Dist_To_Ref"]) as u_cur:
        for row in u_cur:
            oid = row[0]
            row[1] = tf_results.get(oid, 0)
            row[2] = dist_results.get(oid, -1)
            u_cur.updateRow(row)

    arcpy.AddMessage("Analysis complete successfully.")

if __name__ == "__main__":
    main()