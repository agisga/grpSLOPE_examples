output_file_names = (1..20).to_a
output_file_names.map! { |num| "/home/agossman/github/grpSLOPE_examples/EigenPrism/Brzyski_figure_4_with_EigenPrism_#{num}.RData" }

qos     = "normal"
time    = "24:00:00" 
nodes   = 1
ntasks  = 1
mem     = 128000
rscript = "/home/agossman/github/grpSLOPE_examples/Brzyski_figure_4_with_EigenPrism.R"
n_iter  = 10 # number of iterations within each run of rscript (want 200 total iterations)

20.times do |i|
  job_name = "EigenP#{i+1}"

  srun = "#!/bin/bash\n#SBATCH --qos=#{qos}\n#SBATCH --job-name=#{job_name}\n#SBATCH --time=#{time}\n#SBATCH --nodes=#{nodes}\n#SBATCH --ntasks-per-node=#{ntasks}\n#SBATCH --mem=#{mem}\nmodule load R\nmodule load matlab/r2015b\nRscript #{rscript} #{output_file_names[i]} #{n_iter}"

  File.open("Brzyski_figure_4_with_EigenPrism.srun", "w") { |f| f.write(srun) }
  system "sbatch Brzyski_figure_4_with_EigenPrism.srun"
  puts "something went wrong with number #{i+1}" if $?.exitstatus > 0
end


