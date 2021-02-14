#!/usr/bin/env julia

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--opt1"
            help = "an option with an argument"
        "--opt2", "-o"
            help = "another option with an argument"
            arg_type = Int
            default = 0
        "--flag1"
            help = "an option without argument, i.e. a flag"
            action = :store_true
        "arg1"
            help = "a positional argument"
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
end

main()

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
