import FDTD2D
using Base.Test

@time @test FDTD2D.run(Dict("time_steps" => 2, "time_step" => 0.1)) == true
