"""
    set_random(mme::MME,randomStr::AbstractString,vc::Float64, df::Float64))

set variables as iid random effects
"""
function set_random(mme::MME,randomStr::AbstractString, vc::Float64, df::Float64)
    setAsRandom(mme,randomStr, vc, df)
end

"""
    outputMCMCsamples(mme::MME,trmStr::AbstractString...)

Get samples for specific variables.
"""
function outputMCMCsamples(mme::MME,trmStr::AbstractString...)
    for i in trmStr
      outputSamplesFor(mme,i)
    end
end
