#!/usr/bin/env julia

using ArgParse
using QTL

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-i"
            help = "input file name"
            required = true
        "-o"
            help = "output file name"
            required = true
        "-m"
            help = "missing value"
            arg_type = Int
            default = 9
        "--flag1"
            help = "an option without argument, i.e. a flag"
            action = :store_true
    end
    return parse_args(s)
end

function main()
      parsed_args = parse_commandline()

      mygeno=readdlm(parsed_args["i"])
      mymissing=parsed_args["m"]
      missing2mean(mygeno,mymissing)
      writedlm(parsed_args["o"],mygeno," ")
end

main()
