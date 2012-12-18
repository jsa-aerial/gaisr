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
  puts " * embl-to-nc embl-file"
  puts " * entry-file-intersect {'with-seqs'} file1 file2 & files"
  puts " * run-config config-file"
  puts " * check-job {all | jobid+}"
  puts " * rna-taxon-info outfile [rnas | -f rna-file] [taxons | -f taxon-file]"
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
            'embl-to-nc',
            'entry-file-intersect',
            'run-config',
            'check-job',
            'rna-taxon-info'].include?(@cmd))
  puts "Unknown subcmd #{@cmd}"
  puts "use gaisr -h, for useage"
  exit(1)
end


### run-config and check-sto currently local here.  May split these
### out as separate commands as well.


def test_file (f, *msg)
  if File.exists?(f)
    true
  else
    if msg.empty?
      puts "Error: file '#{f}' does not exist"
    else
      puts msg
    end
    exit(1)
  end
end


def remote_cmd (upload_type, infile, *misc)
  server = "http://roz.bc.edu:8080"
  base_url = server + "/mlab/upld"

  if infile.is_a?(Array)
    infile = infile.map do |f| File.new(f) end
  else
    infile = File.new(infile)
  end

  return RestClient.post(base_url,
                         {"upload-type" => upload_type,
                           :subtype => "remote",
                           :user => ENV['LOGNAME'].downcase,
                           :host => ENV["HOST"],
                           :file => infile, ##File.new(infile),
                           :misc => misc,
                           :multipart => true ##NECESSARY when vec of files!
                         },
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

def embl_to_nc_finish(file, info)
  basedir = File.dirname(file)
  filename = File.basename(file)
  filetype = "\.#{filename.split(".").last}"
  puts "Finished EMB to NC conversion for #{file}", "Results in #{basedir}"
  fname = filename.sub(/#{filetype}/, "-NC#{filetype}")
  lines = info[1]
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

  elsif (cmd == "embl-to-nc")
    result = JSON.parse(result.body)
    info = result["info"]
    info_len = info.length
    jstat = info[0]
    if (jstat != "done")
      report_result("", result)
    elsif ((rstat = info[1]) != "good")
      report_result("ERROR:", result)
    else
      embl_to_nc_finish(file, info[2..info_len])
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
    test_file(jobidfile, "Error - you have no jobs")
  else
    jobid = args[0]
    jobidfile = "job"+jobid
    jobidfile = File.join(@gaisrDir, jobidfile)
  end
  if test_file(jobidfile, "No job with id #{jobid} found")
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


def embl_to_nc (args)
  infile = File.expand_path(args[0])
  result = remote_cmd(@cmd, infile)
  result = report_result("", result)
  jobid = result["info"][3]
  jobidfile = create_jobid_file("#{jobid}", [@cmd, infile])
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


def entry_file_intersect(args)
  wsqs = (args[0] == "with-seqs")
  if ((wsqs and args.length < 3) or (not wsqs and args.length < 2))
    puts "entry-file-intersect needs at least two input files"
    exit(1)
  end
  files = if (wsqs) then args[1..args.length] else args end
  result = remote_cmd(@cmd, files, wsqs)
  result = JSON.parse(result.body)
  if (result["stat"] != "success")
    report_result("ERROR:", result)
  else
    info = result["info"][2] # first two are not used remote job status
    basedir = File.dirname(files[0])
    filename, ext = File.basename(files[0]).split(".")
    newfilename = "#{filename}-INTERSECT.#{ext}"
    newfullspec = File.join(basedir, newfilename)
    outFile = File.new(newfullspec, "w")
    sto = (ext == "sto")
    if (!outFile)
      puts "Unable to open file #{newfullspec}"
    else
      info.each do |x|
        if wsqs
          ent, sq = x
          if sto
            outFile.printf("%-35s %-1s\n" % x)
          else
            outFile.printf(">#{ent}\n")
            outFile.printf("#{sq}\n")
          end
        else
          outFile.printf("#{x}\n")
        end
      end
      outFile.close
      puts "Intersection output in #{newfullspec}"
    end
  end
end


def rna_taxon_info (args)
  i = 0
  outname = args[i]

  if (args[i += 1] == "-f")
    if test_file(args[i += 1])
      rnas = IO.readlines(args[i])
    end
  else
    rnas = args[i].split(/ *, */)
  end
  rnas = rnas.map do |r| bits = r.split(/( |\.|-)/); [bits[0], bits[2]] end

  if (args[i += 1] == "-f")
    if test_file(args[i += 1])
      taxons = IO.readlines(args[i])
    end
  else
    taxons = args[i].split(/ *, */)
  end

  time = Time.now
  tempfile = File.join(Dir.tmpdir,
                       ["rnatxinfo", time.hour, time.min, time.sec].join("-"))
  tempFile = File.new(tempfile, "w")
  rnas.each do |(r, v)| tempFile.puts("#{r}-#{v}") end
  tempFile.puts(";;;")
  taxons.each do |t| tempFile.puts(t) end
  tempFile.close
  result = remote_cmd(@cmd, tempfile)
  result = JSON.parse(result.body)
  File.delete(tempfile)

  if (result["stat"] != "success")
    report_result("ERROR:", result)
  else
    outFile = File.new(outname, "w")
    if (!outFile)
      puts "Unable to open file #{outname}"
    else
      info = result["info"][2] # first two are not used remote job status
      outFile.syswrite(info)
      outFile.close
      puts "Taxon group information in #{outname}"
    end
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

  when "embl-to-nc"
    puts "Takes an input file with EMBL entry names and converts it to one"
    puts "with NCBI NC names.  Maintains sequence information.  EMBL-FILE is"
    puts "either a sto or fna file with entries named with EMBL accensions."
    puts "The result file name is the input name with '-NC' appended.  The"
    puts "result file is a corresponding sto or fna with the same sequence"
    puts "(and alignment if sto) with NCBI NC accensions replacing EMBL names."
    puts ""
    puts "Because embl-to-nc may take considerable time (many minutes), it is"
    puts "always run as a job, and so will immediately return the job id for"
    puts "the run.  This may then be checked and results fetched by use of"
    puts "check-job."

  when "entry-file-intersect"
    puts "Intersects two or more sto, aln, fasta, or 'entry' files and writes"
    puts "the result to an output file in the same directory as input and with"
    puts "a filename that is the input name with -INTERSECT appended."
    puts ""
    puts "If 'with-seqs' is given (preceding all input files), the output also"
    puts "has the corresponding sequences for the result entry set.  In this"
    puts "case all input files must be of the same type (all stos, all fnas,"
    puts "etc.) and the file type of the output is the same."
    puts "***NOTE: for sto files, all # lines are removed!"
    puts ""
    puts "If 'with-seqs' is not given, the output is simply an entry file (a"
    puts ".ent file) with entries listed one per line."

  when "rna-taxon-info"
    puts "Perform a 'new rna' taxon grouping information analysis on the given"
    puts "rnas and taxons and place resulting information in the given output"
    puts "file."
    puts ""
    puts "rnas can be given as either a string of comma separated rnas or by"
    puts "using the -f option, can be listed, one per line, in the given file."
    puts "In either case, an rna must be in the format rnaname[.|-| ]version."
    puts "rnaname is a Meyer lab convention name for new rnas, for example,"
    puts "rna_00011.  Version is a positive integer denoting the version of"
    puts "interest.  So, a full example could be rna_00011.2"
    puts ""
    puts "taxons can be given as either a string of comma separated taxon names"
    puts "or by using the -f option, can be listed one per line, in the given"
    puts "file.  The set of default taxons can be denoted by the special name"
    puts "'default-taxons'."
    puts ""
    puts "Examples:"
    puts "gaisr rna-taxon-info rna_11-2-txinfo.txt rna_00011.2 default-taxons"
    puts "gaisr rna-taxon-info taxon-info.txt -f rnas.txt -f taxons.txt"
    puts ""
    puts "Output format:"
    puts "rna-name, version, taxon-name, rna cnt in taxon, total cnt, percent"
    puts "list of genomes (by NC name)"

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

elsif (@cmd == "embl-to-nc")
  embl_to_nc(args)

elsif (@cmd == "entry-file-intersect")
  entry_file_intersect(args)

elsif (@cmd == "rna-taxon-info")
  rna_taxon_info(args)

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
