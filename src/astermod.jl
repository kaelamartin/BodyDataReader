"""
vargout = astermod(aster, prop)

Uses inputs to find requested asteroid shape model properties.
"""
function astermod(aster::Any, prop::String)
    # Asteroid data from Database of Asteroid Models from Inversion Techniques (DAMIT),
    #     http://astro.troja.mff.cuni.cz/projects/asteroids3D/web.php?page=db_export

# Check if asteroid data file exists in current directory
if isfile("asteroids")
    # Open the file
    f = open("asteroids");
    rd = readlines(f);
    close(f)
else
    # Extract updated asteroid data from DAMIT
    ext_data = download("http://astro.troja.mff.cuni.cz/projects/asteroids3D/php/db_export_extended.php", "asteroids");
    # Open the file
    f = open("asteroids");
    rd = readlines(f);
    close(f)
end

# Initialize keys for asteroid properties dictionary, aster_props
k = ["aster_id",
     "aster_num",
     "aster_name",
     "aster_desig",
     "aster_comment",
     "model_id",
     "lambda",
     "beta",
     "period",
     "yorp",
     "jd0",
     "phi0",
     "used_model",
     "model_p1",
     "model_p2",
     "model_p3",
     "model_p4",
     "model_p5",
     "calib_size",
     "equiv_d",
     "equiv_d_err",
     "ref_id",
     "ref_id_obs",
     "therm_i",
     "therm_imin",
     "therm_imax",
     "vis_geo_alb",
     "vis_geo_alb_err",
     "op_crater_ang",
     "a_dens_crater",
     "qual_flag" ,
     "model_ver",
     "model_comment"]

# Initialize model id dictionary
model_id = Dict()

# Iterate through recorded asteroid line number within database
for (i, line_no) in enumerate(rd)

    # Compile properties of asteroid model
    props_parse = collect(m.match for m in eachmatch(r"(\"(.*?)\",)|(\"(.*?)\")|(,)", rd[i]))
    props = [strip(thing, ['\"', '\\', ',']) for thing in props_parse]

    for (i, p) in enumerate(props)
        if p == ""
            props[i] = "Not available."
        else
            continue
        end
    end

    # Create dictionary of asteroid properties
    aster_props = Dict(zip(k,props))

    # Sort by asteroid model
    model_id[props[6]] = aster_props
end

# Initialize model and property list
model_list = []
property_vals = []

for (key, value) in model_id
    # Parse through asteroids file for inputted asteroid
    #print(key)

    if aster == model_id[key]["aster_id"]
        # Checks for input of asteroid id
        aster_property = model_id[key][prop]
        push!(model_list, key)
        push!(property_vals, aster_property)
    elseif aster == model_id[key]["aster_num"]
        # Checks for input of asteroid number
        aster_property = model_id[key][prop]
        push!(model_list, key)
        push!(property_vals, aster_property)
    elseif aster == model_id[key]["aster_name"]
        # Checks for input of asteroid name
        aster_property = model_id[key][prop]
        push!(model_list, key)
        push!(property_vals, aster_property)
    else
        # Throw error if no asteroid exists with inputted method of identification
        #error("ERROR: No asteroid found with inputted identification.")
    end
end

# Output is the model id (key) with the requested property (value)
vargout = Dict(zip(model_list, property_vals))
#print(vargout)
return vargout

    # TODO: add error if fields do not match the database (pairs, etc.)
    # TODO: add capability of requesting multiple properties for a single asteroid
    # TODO: add capability of calling multiple asteroids for a single (or multiple) properties

end
