module FDTD2D
using Plots
export run

"""Calculate parameters needed for simulation and return them as a dictionary."""
function calculate_params(init_params)
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
function init_data(params)
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
function generate_fields_x_min!(d, time, p)
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

"""Update E_x"""
function update_ex!(ex, hz, v)
  for j = 2:size(ex,2)-1
    for i = 2:size(ex,1)-1
      ex[i,j] += v*(hz[i,j] + hz[i-1,j] - hz[i,j-1] - hz[i-1,j-1])
    end
  end
end

"""Update E_y"""
function update_ey!(ey, hz, v)
  for j = 2:size(ey,2)-1
    for i = 2:size(ey,1)-1
      ey[i,j] -= v*(hz[i,j] - hz[i-1,j] + hz[i,j-1] - hz[i-1,j-1])
    end
  end
end

"""Update H_x"""
function update_hx!(hx, ez, v)
  for j = 1:size(hx,2)
    for i = 1:size(hx,1)
      hx[i,j] -= v*(ez[i+1,j+1] + ez[i,j+1] - ez[i+1,j] - ez[i,j])
    end
  end
end

"""Update H_y"""
function update_hy!(hy, ez, v)
  for j = 1:size(hy,2)
    for i = 1:size(hy,1)
      hy[i,j] += v*(ez[i+1,j+1] - ez[i,j+1] + ez[i+1,j] - ez[i,j])
    end
  end
end

"""Make step"""
function make_step!(d, p)
  update_ex!(d["ex"], d["hz"], p["half_cfl"]["y"])
  update_ey!(d["ey"], d["hz"], p["half_cfl"]["x"])
  update_ey!(d["ezx"], d["hy"], -p["half_cfl"]["x"])
  update_ex!(d["ezy"], d["hx"], -p["half_cfl"]["y"])
  d["ez"] = d["ezx"] + d["ezy"]

  update_hx!(d["hx"], d["ez"], p["half_cfl"]["y"])
  update_hy!(d["hy"], d["ez"], p["half_cfl"]["x"])
  update_hy!(d["hzx"], d["ey"], -p["half_cfl"]["x"])
  update_hx!(d["hzy"], d["ex"], -p["half_cfl"]["y"])
  d["hz"] = d["hzx"] + d["hzy"]
end

function output(d, k, p)
  return true
end

"""Run the simulations"""
function run(init_params_dict)
  params_dict = calculate_params(init_params_dict)
  data = init_data(params_dict)
  for k = 1:params_dict["time_steps"]
    time = k*params_dict["time_step"]
    generate_fields_x_min!(data, time, params_dict)
    make_step!(data, params_dict)
    output(data, k, params_dict)
  end
  return data
end

end
