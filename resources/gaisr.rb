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
  puts " * list {['active' | 'done']}"
  puts " * check-sto <stofile>+"
  puts " * get-seqs [sto|aln|fna|ent] {[+|-]prefix [+|-]suffix}"
  puts " * correct-sto-coordinates stofile"
  puts " * embl-to-nc-sto embl-stofile"
  puts " * run-config config-file"
  puts " * check-job {all | jobid+}"
  puts " * rna-taxon-info outfile rna-list taxon-list"
  puts " * name-taxonomy"
  puts " * sto-csv-matchup"
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
            'get-seqs',
            'correct-sto-coordinates',
            'embl-to-nc-sto',
            'run-config',
            'check-job',
            'rna-taxon-info'].include?(@cmd))
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


def get_seqs (args)
  infile = File.expand_path(args[0])
  ld = if (args.length > 1) then args[1] else 0 end
  rd = if (args.length > 2) then args[2] else 0 end
  result = remote_cmd(@cmd, infile, [ld, rd])
  result = JSON.parse(result.body)
  if (result["stat"] == "success")
    info = result["info"]
    basedir = File.dirname(infile)
    filename = File.basename(infile)
    newfilename = "new-#{filename}"
    newfullspec = File.join(basedir, newfilename)
    outFile = File.new(newfullspec, "w")
    if outFile
      outFile.syswrite(info[info.length-1])
      outFile.close
      puts "New/Adjusted entries/seqs in #{newfullspec}"
    else
      puts "Unable to open file #{newfilename}"
    end
  else
    report_result("ERROR:", result)
  end
end




def display_help (cmd)
  puts ""
  case cmd
  when "list"
    puts "Lists your current ('active') jobs or your set of completed jobs."
    puts "With 'active' given, lists active/current running jobs."
    puts "With 'done' given, lists the completed jobs."
    puts "Default, when no argument is given, lists active jobs."

  when "check-sto"
    puts "Runs the gaisr sto file integrity checker.  Checks for a legally"
    puts "constructed and formatted file.  If the file is 'good', returns"
    puts "that as a message.  If there are problems, the problems are listed"
    puts "along with the genome names (N numbers) of the offending lines."

  when "get-seqs"
    puts "Takes a sto, aln, fasta, or 'entry' file and optional integers"
    puts "indicating how much to add or subtract to the start and end of"
    puts "sequences. Negative values always trim while positive values always"
    puts "add elements (nucleotides) to the sequence.\n"
    puts ""
    puts "Entry files are text files which have one entry per line, where an"
    puts "entry is <name/coordinates/strand>. 'name' is a genome name (only"
    puts "NC_* names currently).  'coordinates' are start-end, where start"
    puts "and end are integers, with end > start.  Strand is 1 or -1, and"
    puts "indicates the strand. For example: NC_013525/2 1770792-1770817/1\n"
    puts ""
    puts "get-seqs Supports two basic ways of getting or adjusting sequences.\n"
    puts ""
    puts "* For entry files and fasta files, just takes the entry coordinates,"
    puts "  adjusts them by any given prefix/suffix, and returns the sequence"
    puts "  for the updated coordinates (for entry name and strand) in fasta"
    puts "  format."
    puts ""
    puts "* For alignment files (sto and clustalw aln) the coordinates of"
    puts "  entries along with the given prefix/suffix are used to 'adjust'"
    puts "  the aligned (gapped) corresponding sequences.  For this to work"
    puts "  the coordinates _must_ correctly match exactly the degapped"
    puts "  version of the provided sequence."
    puts ""
    puts "  To ensure this is the case, you can use the correct-sto-coordinates"
    puts "  gaisr command to obtain a version of a sto file with each entry's"
    puts "  coordinates corrected to ensure a match with its given sequence."
    puts ""
    puts "  Gaps are treated as characters and counted in terms of trimming the"
    puts "  sequences (removal of characters from either end), but are not"
    puts "  counted in adjusting the corresponding entry's coordinates.  So,"
    puts "  the count of non gap characters in the trimmed subsequence is used"
    puts "  to compute the new coordinates that match the trimmed sequence."

  when "correct-sto-coordinates"
    puts "Takes a sto file, and for each entry/sequence pair, determines if the"
    puts "degapped sequence occurs _exactly_ in the genome given in the entry."
    puts "An entry is as defined for get-seqs.  Because correct-sto-coordinates"
    puts "may take considerable time (many minutes), it is always run as a job,"
    puts "and so will immediately return the job id for the run.  This may then"
    puts "be checked and results fetched by use of check-job.  Returns up to"
    puts "3 new files:"
    puts ""
    puts "* A new version of the original input sto file with the same name as"
    puts "  the original, with the original renamed with a trailing -old."
    puts ""
    puts "* A differences file, named off the original sto file name with a"
    puts "  trailing -diffs.  This contains an entry mapping from any entries"
    puts "  in the original to corresponding new versions as determined by the"
    puts "  true location and strand of the sequence in the genome."
    puts ""
    puts "* A 'bad' file, named off the original sto file name with a trailing"
    puts "  -bad.  This contains entries whose sequences could not be found in"
    puts "  the entry's given genome.  Each 'bad' entry has the original entry"
    puts "  (name/coordinates/strand), the original gapped sequence given for"
    puts "  the entry, the degapped version of this sequence, and the actual"
    puts "  sequence in the genome at the given coordinates."
    puts ""
    puts "If no differences or bad entries were found, the diffs and bad files"
    puts "will not exist, however there is always a new version file with the"
    puts "original renamed to old."

  when "run-config"
    puts "Takes a job configuration (job config or simply config) file and"
    puts "submits it for running as a background job.  The config file is"
    puts "documented with config -h, and encompasses many different tools and"
    puts "options (such as cmbuild, cmcalibrate, cmsearch, csv-gen, et.al.)"
    puts "These are typically quite long running (hours to even days), and"
    puts "run-config returns immediately the job id for the job.  This can be"
    puts "checked, and results obtained with check-job."

  when "check-job"
    puts "Checks on a job that was submitted for running in Gaisr in the"
    puts "background.  The job may have been explicitly submitted by run-config"
    puts "or implicitly via another command such as correct-sto-coordinates."
    puts "Takes either the keyword 'all' or a list of one or more jobids (as"
    puts "returned from a job submitting command)"
    puts ""
    puts "For each requested job (or all current jobs if 'all' is given),"
    puts "check-job inquires of Gaisr the current status of the job.  If it is"
    puts "finished, it obtains the results of the job, otherwise reports the"
    puts "current running status of the job."
  end
  puts ""
end




if ((args.length > 0) and (["-h", "--help"].include?(args[0])))
    display_help(@cmd)

elsif (@cmd == "list")
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

elsif (@cmd == "get-seqs")
  if (args.length == 0)
    puts "#{@cmd} requires at least a file."
    exit(1)
  end
  get_seqs(args)

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
