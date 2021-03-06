include("../src/FDTD2D.jl")
import .FDTD2D
using Plots
using Test

# Test FDTD2D.calculate_params
function laser_pulse_gauss(x_min, width, duration)
  function f(t, x, y)
    return exp(-(y/width)^2)*exp(-((t - (x - x_min) - duration)/duration)^2)*sin(2.0*pi*(t - x))
  end
  return f
end

params_correct = Dict(
  "x_bounds" => (0.0, 12.0),
  "y_bounds" => (-4.0, 4.0),
  "matrix_size" => Dict("x" => 96, "y" => 64),
  "time_steps" => 81,

  "output" => Dict(
    "iteration_pass" => 10,
    "directory_name" => "test/output"
  )
)
params_correct["laser_pulse_y_shape"] = laser_pulse_gauss(params_correct["x_bounds"][1], 1., 2.)
params_correct["laser_pulse_z_shape"] = laser_pulse_gauss(params_correct["x_bounds"][1], 1., 2.)

params_calculated = FDTD2D.calculate_params(params_correct)

@testset "Calculation of parameters" begin
    @test params_calculated["box_size"]["x"] == 12.0
    @test params_calculated["box_size"]["y"] == 8.0
    @test params_calculated["space_step"]["x"] == 0.125
    @test params_calculated["space_step"]["y"] == 0.125
    @test params_calculated["time_step"] == sqrt(2.0)/16.0
    @test params_calculated["half_cfl"]["x"] == 0.5*sqrt(0.5)
    @test params_calculated["half_cfl"]["y"] == 0.5*sqrt(0.5)
end

# Test FDTD2D.init_data
data_empty = Dict()
for name in ("ex","ey","ezx","ezy","ez")
  data_empty[name] = zeros(Float64, params_calculated["matrix_size"]["x"], params_calculated["matrix_size"]["y"])
end
for name in ("hx","hy","hzx","hzy","hz")
  data_empty[name] = zeros(Float64, params_calculated["matrix_size"]["x"]-1, params_calculated["matrix_size"]["y"]-1)
end

data_test = FDTD2D.init_data(params_calculated)
@test isequal(data_test,data_empty)

@time d = FDTD2D.run(params_correct)

heatmap(d["ez"])
