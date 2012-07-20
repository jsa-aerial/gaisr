#!/usr/bin/env ruby

require 'rubygems'
require 'rest_client'
require 'json'
require 'fileutils'
require 'tmpdir'

#ARGV.each do|x|
#  puts "Argument: #{x}"
#end

@gaisrDir = File.expand_path("~/.gaisr")
@doneDir = File.join(@gaisrDir, "Finished")

if (!File.directory?(@gaisrDir))
  Dir.mkdir(@gaisrDir, 0775)
end

if (!File.directory?(@doneDir))
  Dir.mkdir(@doneDir, 0775)
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
  puts "\nUsage:\n"
  puts "gaisr subcmd args"
  puts ""
  puts "subcmd is one of\n\n"
  puts " * list <'active'> | <'done'>"
  puts " * check-sto <stofile>+"
  puts " * add-sequence <stofile> [+|-]prefix [+|-]suffix"
  puts " * correct-sto-coordinates <stofile>"
  puts " * run-config <config-file>"
  puts " * check-job {all | <jobid>+}"
  puts " * name-taxonomy\n * sto-csv-matchup"
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
elsif (not ['list',
            'check-sto',
            'add-sequence',
            'correct-sto-coordinates',
            'run-config',
            'check-job'].include?(@cmd))
  puts "Unknown subcmd #{@cmd}"
  puts "use gaisr -h, for useage"
  exit(1)
end


### run-config and check-sto currently local here.  May split these
### out as separate commands as well.


def remote_cmd (upload_type, infile, *misc)
  server = "http://roz.bc.edu:8080"
  base_url = server + "/mlab/upld"
  return RestClient.post(base_url,
                         {"upload-type" => upload_type,
                           :subtype => "remote",
                           :user => ENV['LOGNAME'].downcase,
                           :host => ENV["HOST"],
                           :file => File.new(infile),
                           :misc => misc},
                         {:cookies => {:user => ENV['LOGNAME'].downcase}})
end


def report_result (header, result)
  result = if result.is_a?(String) then JSON.parse(result.body) else result end
  if (result["stat"] != "success")
    puts "Error - ", result["info"]
  else
    puts ""
    if (header != "")
      puts header
    end
    info = result["info"]
    info_len = info.length
    jstat = info[0]
    ##puts "JSTAT = #{jstat}"
    ##puts "info[1..] = #{info[1..info_len]}"
    msg_bits = info[2..info_len]
    if ((jstat == "done") and ((rstat = info[1]) != "good"))
      puts "ERROR:"
    end
    puts msg_bits
  end
  result
end


def create_jobid_file (jobid, content)
  jobidfile = @gaisrDir + "/job"+jobid
  jobidFile = File.new(jobidfile, "w")
  jobidFile.syswrite("#{jobid}\n")
  if (!content.empty?)
    jobidFile.syswrite(content.join("\n"))
  end
  jobidFile.close
  FileUtils.cp(jobidfile, @gaisrDir + "/last-job.txt")
  return jobidfile
end


def correct_sto_finish (file, info)
  ## it's amazing how verbose and convoluted this all is.  You
  ## could easily do this entire function in three lines of Clojure.
  ## Ruby sucks!!!
  basedir = File.dirname(file)
  filename = File.basename(file)
  puts "Finished sto correction for #{file}", "Results in #{basedir}"

  contents = info[1..info.length]
  suffixes = ["-new.sto", "-diffs.txt", "-bad.txt"]
  names = suffixes.map do |x| filename.sub(/\.sto/, x) end
  both = [names, contents]
  both.transpose.each do |fname, lines|
      if (!lines.empty?)
        outFile = File.new(File.join(basedir, fname), "w")
        if outFile
          outFile.syswrite(lines)
          outFile.close
        else
          puts "Unable to open file #{fname}"
        end
      end
    end
  if File.exist?(File.join(basedir, names[0]))
    File.rename(file, File.join(basedir, filename.sub(/\.sto/, "-old.sto")))
    File.rename(File.join(basedir, names[0]), file)
  end
end


def check_finish (jobidfile, data, result)
  jobid, cmd, file = data
  #puts "Data ", data, jobid, cmd, file
  if (cmd == "run-config")
    result = report_result("", result)
    info = result["info"]
    jstat = info[0]
  elsif (cmd == "correct-sto-coordinates")
    result = JSON.parse(result.body)
    info = result["info"]
    info_len = info.length
    jstat = info[0]
    if (jstat != "done")
      report_result("", result)
    elsif ((rstat = info[1]) != "good")
      report_result("ERROR:", result)
    else
      correct_sto_finish(file, info[2..info_len])
    end
  end
  if ((jstat == "done") and File.exist?(jobidfile))
    filename = File.basename(jobidfile)
    File.rename(jobidfile, File.join(@doneDir, filename))
  end
end


def check_job (args)
  if (args == [])
    jobidfile = File.join(@gaisrDir, "last-job.txt")
    if (!File.exists?(jobidfile))
      puts "Error - you have no jobs"
      exit(1)
    end
  else
    jobid = args[0]
    jobidfile = "job"+jobid
    jobidfile = File.join(@gaisrDir, jobidfile)
  end
  if File.exists?(jobidfile)
    data = IO.readlines(jobidfile).map do |x| x.chomp end
    jobid = data[0]
    tempfile = File.join(Dir.tmpdir, jobid)
    tempFile = File.new(tempfile, "w")
    tempFile.puts(jobid)
    tempFile.close
    result = remote_cmd(@cmd, tempfile)
    File.delete(tempfile)
    ## Make sure we pass in the actual jobidfile (not default last-job.txt!!)
    jobidfile = File.join(@gaisrDir, "job"+jobid)
    check_finish(jobidfile, data, result)
  else
    puts "No job with id #{jobid} found"
  end
end

def check_jobs(jobids)
  jobids.each do |x|
    check_job([x])
  end
end


def correct_sto (args)
  args.each do |infile|
    infile = File.expand_path(infile)
    result = remote_cmd(@cmd, infile)
    result = report_result("", result)
    jobid = result["info"][3]
    jobidfile = create_jobid_file("#{jobid}", [@cmd, infile])
  end
end


def run_job (config_file)
  result = remote_cmd(@cmd, config_file)
  result = report_result("", result)
  jobid = result["info"][3]
  jobidfile = create_jobid_file("#{jobid}", [@cmd, config_file])
end


def list_jobs (state)
  dir = if (state == "active") then @gaisrDir else @doneDir end
  curjobs = Dir.glob(File.join(dir, "job*"))
  jstate = if (state == "active") then "running" else "finished" end
  if (curjobs.length == 0)
    puts "You have no #{jstate} jobs..."
  else
    curjobs.each do |jobfile|
      jobid, cmd, file = IO.readlines(jobfile).map do |x| x.chomp end
      puts "Job #{jobid}, #{jstate} command #{cmd} on #{file}"
    end
  end
end


if (@cmd == "list")
  if ((args.length == 0) or (["active", "running"].include?(args[0])))
    list_jobs("active")
  elsif (["done", "finished"].include?(args[0]))
    list_jobs("done")
  else
    puts "Unknown list request '#{args[0]}'"
  end

elsif (@cmd == "check-sto")
  args.each do |infile|
    infile = File.expand_path(infile)
    result = remote_cmd(@cmd, infile)
    report_result(infile, result)
  end

elsif (@cmd == "correct-sto-coordinates")
  correct_sto(args)

elsif (@cmd == "check-job")
  if ((args.length != 0) and (args[0] == "all"))
    curjobs = Dir.glob(File.join(@gaisrDir, "job*"))
    if (curjobs.length == 0)
      puts "You have no active jobs..."
    else
      jobids = curjobs.map do |x| IO.readlines(x)[0].chomp end
      check_jobs(jobids)
    end
  elsif (args.length > 1)
    check_jobs(args)
  else
    check_job(args)
  end

elsif (args.length > 1)
  puts "#{@cmd} has only single config file parameter, not #{args.length}"
  exit(1)
else
  config_file = File.expand_path(args[0])
  run_job(config_file)
end
