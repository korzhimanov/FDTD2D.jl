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

  params["y1"] = range(params["y_bounds"][1], length=params["matrix_size"]["y"], step=params["space_step"]["y"]) |> collect
  params["y2"] = range(params["y_bounds"][1]+0.5*params["space_step"]["y"], length=params["matrix_size"]["y"]-1, step=params["space_step"]["y"]) |> collect

  return params
end

"""Initialize data"""
function init_data(params::Dict)
  data = Dict{String, Array{Float64,2}}()
  for name in ("ex","ey","ezx","ezy","ez")
    data[name] = zeros(Float64, params["matrix_size"]["x"], params["matrix_size"]["y"])
  end
  for name in ("hx","hy","hzx","hzy","hz")
    data[name] = zeros(Float64, params["matrix_size"]["x"]-1, params["matrix_size"]["y"]-1)
  end
  return data
end

"""Generate fields at left boundary (x_min = p["x_bounds"][0]) at time moment 'time'"""
function generate_fields_x_min!(d::Dict, time::Float64, p::Dict)
  for j = 1:p["matrix_size"]["y"]
    d["ey"][2,j] -= p["laser_pulse_y_shape"](time, p["x_bounds"][1] + p["space_step"]["x"], p["y1"][j])
    d["ezx"][2,j] -= p["laser_pulse_z_shape"](time, p["x_bounds"][1] + p["space_step"]["x"], p["y1"][j])
  end

  for j = 1:p["matrix_size"]["y"]-1
    d["hzx"][2,j] -= p["laser_pulse_y_shape"](time + 0.5*p["time_step"], p["x_bounds"][1] + 1.5*p["space_step"]["x"], p["y2"][j])
    d["hy"][2,j] -= p["laser_pulse_z_shape"](time + 0.5*p["time_step"], p["x_bounds"][1] + 1.5*p["space_step"]["x"], p["y2"][j])
  end

  d["ez"][2,:] = d["ezx"][2,:] + d["ezy"][2,:]
  d["hz"][2,:] = d["hzx"][2,:] + d["hzy"][2,:]
end

"""Make step"""
function make_step!(d::Dict, p::Dict)
  for j = 2:p["matrix_size"]["y"]-1
    for i = 2:p["matrix_size"]["x"]-1
      d["ex"][i,j] += p["half_cfl"]["y"]*(d["hz"][i,j] + d["hz"][i-1,j] - d["hz"][i,j-1] - d["hz"][i-1,j-1])
      d["ey"][i,j] -= p["half_cfl"]["x"]*(d["hz"][i,j] + d["hz"][i,j-1] - d["hz"][i-1,j] - d["hz"][i-1,j-1])
      d["ezx"][i,j] += p["half_cfl"]["x"]*(d["hy"][i,j] + d["hy"][i,j-1] - d["hy"][i-1,j] - d["hy"][i-1,j-1])
      d["ezy"][i,j] -= p["half_cfl"]["y"]*(d["hx"][i,j] + d["hx"][i-1,j] - d["hx"][i,j-1] - d["hx"][i-1,j-1])
      d["ez"][i,j] = d["ezx"][i,j] + d["ezy"][i,j]
    end
  end

  for j = 1:p["matrix_size"]["y"]-1
    for i = 1:p["matrix_size"]["x"]-1
      d["hx"][i,j] -= p["half_cfl"]["y"]*(d["ez"][i+1,j+1] + d["ez"][i,j+1] - d["ez"][i+1,j] - d["ez"][i,j])
      d["hy"][i,j] += p["half_cfl"]["x"]*(d["ez"][i+1,j+1] + d["ez"][i+1,j] - d["ez"][i,j+1] - d["ez"][i,j])
      d["hzx"][i,j] -= p["half_cfl"]["x"]*(d["ey"][i+1,j+1] + d["ey"][i+1,j] - d["ey"][i,j+1] - d["ey"][i,j])
      d["hzy"][i,j] += p["half_cfl"]["y"]*(d["ex"][i+1,j+1] + d["ex"][i,j+1] - d["ex"][i+1,j] - d["ex"][i,j])
      d["hz"][i,j] = d["hzx"][i,j] + d["hzy"][i,j]
    end
  end
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
