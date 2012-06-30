#!/usr/bin/env ruby

require 'rubygems'
require 'rest_client'
require 'json'

server = "http://roz.bc.edu:8080"
base_url = server + "/mlab/upld"

#puts base_url

#ARGV.each do|x|
#  puts "Argument: #{x}"
#end

if (ARGV.length() == 0)
  puts "use -h for usage"
  exit(1)
end

arg1= ARGV[0]
if (arg1 == "-h" or arg1 == "--help")
  puts "name-taxonomy infilespec outfilespec"
  puts "Generates a taxonomy file from the genome entries in infilespec"
  puts "infilespec can be a sto, fasta, or txt with genome name per line"
  puts ""
  puts "Reports 'success' or if failure, the exception text of failure"
  exit(0)
elsif (ARGV.length() < 2)
  puts "Need to supply 2 arguments: infile outfile"
  puts "use name-taxonomy -h, for useage"
  exit(1)
elsif (ENV["HOST"] != "roz")
  ## puts "Using remote access to mlab..."
end

infile = File.expand_path(ARGV[0])
otfile = File.expand_path(ARGV[1])

result = RestClient.post(base_url,
                         {"upload-type" => 'name2tax',
                           :subtype => "remote",
                           :user => ENV['LOGNAME'],
                           :file => File.new(infile)},
                         {:cookies => {:user => ENV['LOGNAME']}})
result = JSON.parse(result.body)

if (result["stat"] != "success")
  puts "Error - #{result['stat']}"
else
  outFile = File.new(otfile, "w")
  if outFile
    outFile.syswrite(result["info"])
    outFile.close
    puts "success"
  else
    puts "Unable to open file #{otfile}"
  end
end

