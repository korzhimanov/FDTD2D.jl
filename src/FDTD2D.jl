module FDTD2D
export run

function calculate_params(init_params::Dict)
  return init_params
end

function init_data(params::Dict)
  return Dict{String, Number}()
end

"""Generate fields at left boundary (x_min = p['x_bounds'][0]) at time moment 'time'"""
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
