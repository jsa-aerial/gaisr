#!/usr/bin/env ruby

require 'rubygems'
require 'rest_client'
require 'json'

#ARGV.each do|x|
#  puts "Argument: #{x}"
#end

if (ARGV.length() == 0)
  puts "Use -h or --help for usage"
  exit(1)
end

sub_cmds = ["name-taxonomy",
            "sto-csv-matchup"]

cmd= ARGV[0]
args = ARGV[1,ARGV.length]

if (cmd == "-h" or cmd == "--help")
  puts "Usage:\n"
  puts "gaisr subcmd args"
  puts ""
  puts "subcmd is one of\n * name-taxonomy\n * sto-csv-matchup\n * run-config"
  puts " * check-sto\n * correct-sto-coordinates\n"
  puts ""
  puts "args is the set of arguments appropriate to subcmd (including -h)"
  puts ""
  puts "Reports the results of subcmd's execution.  Typically 'success'"
  puts "or if failure, the exception text of failure"
  exit(0)
end

if (sub_cmds.include?(cmd))
  argstr = args.join(" ")
  puts `#{cmd} #{argstr}` # execute cmd, results go to term
  exit(0)
elsif (not ['check-sto', 'run-config'].include?(cmd))
  puts "Unknown subcmd #{cmd}"
  puts "use gaisr -h, for useage"
  exit(1)
end


### run-config and check-sto currently local here.  May split these
### out as separate commands as well.


def remote_cmd (upload_type, infile)
  server = "http://roz.bc.edu:8080"
  base_url = server + "/mlab/upld"
  return RestClient.post(base_url,
                         {"upload-type" => upload_type,
                           :subtype => "remote",
                           :user => ENV['LOGNAME'],
                           :host => ENV["HOST"],
                           :file => File.new(infile)},
                         {:cookies => {:user => ENV['LOGNAME']}})
end

if (cmd == "check-sto")
  args.each do |infile|
    ##infile = File.expand_path(args[0])
    result = remote_cmd(cmd, infile)
    result = JSON.parse(result.body)
    if (result["stat"] != "success")
      puts "Error - #{result['stat']}"
    else
      puts ""
      puts infile
      puts result["info"]
    end
  end
elsif (args.length > 1)
  puts "#{cmd} has only single config file parameter, not #{args.length}"
  exit(1)
else
  result = remote_cmd(cmd, infile)
  result = JSON.parse(result.body)
    if (result["stat"] != "success")
      puts "Error - #{result['stat']}"
    else
      puts "NYI"
    end
end


