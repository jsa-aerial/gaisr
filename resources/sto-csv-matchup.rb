#!/usr/bin/env ruby

require 'rubygems'
require 'rest_client'
require 'json'

server = "http://roz.bc.edu:8080"

#puts base_url

#ARGV.each do|x|
#  puts "Argument: #{x}"
#end

if (ARGV.length() == 0)
  puts "Use -h or --help for usage"
  exit(1)
end

type = "local"
base_url = server + "/mlab/action"
arg1= ARGV[0]

if (arg1 == "-h" or arg1 == "--help")
  puts "sto-csv-matchup stofilespec csvfilespec outfilespec"
  puts ""
  puts "stofilespec can be a sto, fasta, or txt with genome name per line"
  puts ""
  puts "Reports 'success' or if failure, the exception text of failure"
  exit(0)
elsif (ARGV.length() < 3)
  puts "Need to supply 3 arguments: infile csvfile outfile"
  puts "use sto-csv-matchup -h, for useage"
  exit(1)
elsif (ENV["HOST"] != "roz")
  type = "remote"
  base_url = server + "/mlab/upld"
end

sto = File.expand_path(ARGV[0])
csv = File.expand_path(ARGV[1])
out = File.expand_path(ARGV[2])


def local_sto_csv ()
  result = RestClient.post(base_url,
                           {:user => ENV['LOGNAME'],
                             :act => "stocsv-match",
                             :subtype => type,
                             :sto => sto,
                             :csv => csv,
                             :out => out},
                           {:cookies => {:user => ENV['LOGNAME']}})
  result = JSON.parse(result.body)
  if (result["stat"] != "success")
    return "Error - #{result['stat']}"
  else
    return "success"
  end
end


def remote_sto_csv ()
  result = RestClient.post(base_url,
                           {"upload-type" => "stocsv-match",
                             :subtype => type,
                             :user => ENV['LOGNAME'],
                             :file => File.new(sto),
                             :csv => File.new(csv),
                             :out => out},
                           {:cookies => {:user => ENV['LOGNAME']}})
  result = JSON.parse(result.body)
  if (result["stat"] != "success")
    return "Error - #{result['stat']}"
  else
    outFile = File.new(out, "w")
    if outFile
      outFile.syswrite(result["info"])
      outFile.close
      return "success"
    else
      return "Unable to open file #{out}"
    end
  end
end


if (type == "remote")
  result = remote_sto_csv()
else
  result = local_sto_csv()
end

puts result
