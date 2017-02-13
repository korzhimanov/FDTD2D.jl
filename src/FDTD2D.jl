module FDTD2D
export run

"""Calculate parameters needed for simulation and return them as a dictionary."""
function calculate_params(init_params::Dict)
  params = deepcopy(init_params)

  params["box_size"] = Dict{String,Float64}()
  params["box_size"]["x"] = params["x_bounds"][2] - params["x_bounds"][1]
  params["box_size"]["y"] = params["y_bounds"][2] - params["y_bounds"][1]

  params["space_step"] = Dict{String,Float64}()
  params["space_step"]["x"] = params["box_size"]["x"]/params["matrix_size"]["x"]
  params["space_step"]["y"] = params["box_size"]["y"]/params["matrix_size"]["y"]

  params["time_step"] = 0.5*sqrt(params["space_step"]["x"]^2 + params["space_step"]["y"]^2)

  params["half_cfl"] = Dict{String,Float64}()
  params["half_cfl"]["x"] = 0.5*params["time_step"]/params["space_step"]["x"]
  params["half_cfl"]["y"] = 0.5*params["time_step"]/params["space_step"]["y"]

  params["y1"] = [0.0:params["matrix_size"]["y"]-1.0] * params["space_step"]["y"] + params["y_bounds"][1]
  params["y2"] = [0.5:params["matrix_size"]["y"]-1.5] * params["space_step"]["y"] + params["y_bounds"][1]

  return params
end

function init_data(params::Dict)
  return Dict{String, Number}()
end

"""Generate fields at left boundary (x_min = p["x_bounds"][0]) at time moment 'time'"""
function generate_fields_x_min!(d::Dict, time::Float64, p::Dict)
  return true
end

function make_step!(d::Dict, p::Dict)
  return true
end

function output(d::Dict, k::Int, p::Dict)
  return true
end

"""Run the simulations"""
function run(init_params_dict::Dict)
  params_dict = calculate_params(init_params_dict)
  data = init_data(params_dict)
  for k = 1:params_dict["time_steps"]
    time = k*params_dict["time_step"]
    generate_fields_x_min!(data, time, params_dict)
    make_step!(data, params_dict)
    output(data, k, params_dict)
  end
  return true
end

end
