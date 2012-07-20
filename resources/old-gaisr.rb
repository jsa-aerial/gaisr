#!/usr/bin/env ruby

require 'rubygems'
require 'rest_client'
require 'json'
require 'fileutils'

#ARGV.each do|x|
#  puts "Argument: #{x}"
#end

@gaisrDir = File.expand_path("~/.gaisr")

if (!File.directory?(@gaisrDir))
  Dir.mkdir(@gaisrDir)
end

if (ARGV.length() == 0)
  puts "Use -h or --help for usage"
  exit(1)
end

sub_cmds = ["name-taxonomy",
            "sto-csv-matchup"]

@cmd= ARGV[0]
args = ARGV[1,ARGV.length]

if (@cmd == "-h" or @cmd == "--help")
  puts "Usage:\n"
  puts "gaisr subcmd args"
  puts ""
  puts "subcmd is one of\n * name-taxonomy\n * sto-csv-matchup\n * check-sto"
  puts " * correct-sto-coordinates\n * run-config\n * check-job"
  puts ""
  puts "args is the set of arguments appropriate to subcmd (including -h)"
  puts ""
  puts "Reports the results of subcmd's execution.  Typically 'success'"
  puts "or if failure, the exception text of failure"
  exit(0)
end

if (sub_cmds.include?(@cmd))
  argstr = args.join(" ")
  puts `#{@cmd} #{argstr}` # execute cmd, results go to term
  exit(0)
elsif (not ['check-sto', 'correct-sto-coordinates',
            'run-config', 'check-job'].include?(@cmd))
  puts "Unknown subcmd #{@cmd}"
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
                           :user => ENV['LOGNAME'].downcase,
                           :host => ENV["HOST"],
                           :file => File.new(infile)},
                         {:cookies => {:user => ENV['LOGNAME'].downcase}})
end

def report_result (header, result)
  result = JSON.parse(result.body)
  if (result["stat"] != "success")
    puts "Error - #{result['stat']}"
  else
    puts ""
    puts header
    puts result["info"]
  end
  result
end


def check_job (args)
  if (args == [])
    jobidfile = @gaisrDir + "/last-job.txt"
    if (!File.exists?(jobidfile))
      puts "Error - you have no jobs"
      exit(1)
    end
  else
    jobid = args[0]
    jobidfile = "job"+jobid
    jobidfile = @gaisrDir + "/" + jobidfile
  end
  if File.exists?(jobidfile)
    result = remote_cmd(@cmd, jobidfile)
    report_result("Status :", result)
  else
    puts "No job with id #{jobid} found"
  end
end


def run_job (config_file)
  result = remote_cmd(@cmd, config_file)
  result = report_result("", result)
  jobid = result["info"][1]
  jobidfile = @gaisrDir + "/job"+jobid
  jobidFile = File.new(jobidfile, "w")
  jobidFile.syswrite("#{jobid}\n")
  jobidFile.close
  FileUtils.cp(jobidfile, @gaisrDir + "/last-job.txt")
end


if ((@cmd == "check-sto") or (@cmd == "correct-sto-coordinates"))
  args.each do |infile|
    infile = File.expand_path(infile)
    result = remote_cmd(@cmd, infile)
    report_result(infile, result)
  end

elsif (@cmd == "check-job")
  check_job(args)

elsif (args.length > 1)
  puts "#{@cmd} has only single config file parameter, not #{args.length}"
  exit(1)
else
  config_file = File.expand_path(args[0])
  run_job(config_file)
end


