using JWAS
using Base.Test

# write your own tests here
@test 1 == 1
@test ?JWAS.MT.runMCMC()
@test ?JWAS.ST.set_random()
