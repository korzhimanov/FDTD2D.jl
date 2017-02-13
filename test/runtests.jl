import FDTD2D
using Base.Test

function laser_pulse_gauss(x_min::Float64, width::Float64, duration::Float64)
  function f(t::Float64, x::Float64, y::Float64)
    return exp(-(y/width)^2)*exp(-((t - (x - x_min))/duration-2.0)^2)*sin(2.*pi*(t - x))
  end
  return f
end

params_correct = Dict(
  "x_bounds" => (0.0, 8.0),
  "y_bounds" => (-4.0, 4.0),
  "matrix_size" => Dict("x" => 64, "y" => 64),
  "time_steps" => 81,

  "output" => Dict(
    "iteration_pass" => 10,
    "directory_name" => "test/output"
  )
)
params_correct["laser_pulse_y_shape"] = laser_pulse_gauss(params_correct["x_bounds"][1], 1., 2.)
params_correct["laser_pulse_z_shape"] = laser_pulse_gauss(params_correct["x_bounds"][1], 1., 2.)

params_calculated = FDTD2D.calculate_params(params_correct)

@test params_calculated["box_size"]["x"] == 8.0
@test params_calculated["box_size"]["y"] == 8.0
@test params_calculated["space_step"]["x"] == 0.125
@test params_calculated["space_step"]["y"] == 0.125
@test params_calculated["time_step"] == sqrt(2.0)/16.0
@test params_calculated["half_cfl"]["x"] == 0.5*sqrt(0.5)
@test params_calculated["half_cfl"]["y"] == 0.5*sqrt(0.5)

@time @test FDTD2D.run(params_correct) == true
