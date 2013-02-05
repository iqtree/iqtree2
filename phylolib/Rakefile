# recompile and run raxml with rake...

# define what you wanna play with...
VALID_SETS = {
  '7aa' => {:treefile => "",       :dataset => "7.aa.singlegene.binary"}, # AA support?
  '7' => {:treefile => "",       :dataset => "7.dna.singlegene.binary"},
  '20' => {:treefile => "intree20",       :dataset => "20.dna.singlegene.binary"},
  '20part' => {:treefile => "intree20",   :dataset => "20.dna.binary"},   # partition support?
  '49' =>{:treefile => "intree49",        :dataset => "49.binary"},
  '1288' =>{:treefile => "intree1288",    :dataset => "1288.binary"},
} 
VERSIONS = {
  "std" => {:bin => "raxmlLight", :makefile => "Makefile.SSE3.gcc"},
  "bayes" => {:bin => "raxmlBayes", :makefile => "Makefile.BAYES.SSE3.gcc"}
}
STARTING_TREE_OPTS = %w(treefile random parsimony)
VALID_MODELS = %w(GAMMA PSR)
DATADIR="testdata"

# helpers
def print_usage 
    puts "rake run[dataset,version,model,starting_tree,generations]"
    puts "Valid datasets are #{VALID_SETS.keys.join(",")}, expected in folder #{DATADIR}"
    puts "Valid versions are #{VERSIONS.keys.join(",")}"
    puts "Valid models are #{VALID_MODELS.join(",")}"
    puts "Valid starting_trees are #{STARTING_TREE_OPTS.join(",")}"
    puts "generations is only relevant for the bayesian case"
end

def inform_and_delete_file(filename, silent = false)
  if File.exists?(filename) then
    puts "Deleting #{filename}" unless silent
    File.delete(filename)
  else
    puts "?? #{filename} not found"
  end
end

def inform_and_delete_extension_files(extension, silent = false)
  Dir.entries(Dir.pwd).select{|f| File.extname(f) == extension}.each do |f| 
    inform_and_delete_file(f, silent)
  end
end

def run_raxml(version, set, model, starting_tree, generations)
  begin
    raise "invalid version #{version}" unless VERSIONS.keys.include?(version)
    raise "invalid dataset #{set}" unless VALID_SETS.keys.include?(set)
    treefile = VALID_SETS[set][:treefile] 
    dataset =  VALID_SETS[set][:dataset]
    raxmlbin = VERSIONS[version][:bin]
    name = [dataset, version, starting_tree, model].join("_")
    if File.exists?(raxmlbin) 
      tree = File.join DATADIR, treefile
      data = File.join DATADIR, dataset
      call = "./#{raxmlbin} -s #{data} -m #{model}"
      call += " -n #{name}" 
      case starting_tree
        when "treefile" then call += " -t #{tree}"
        when "parsimony" then call += " -p "
        else puts "using random starting tree!"
      end
      if version == "bayes"
        call += " -b -g #{generations}"  
      else
        call = "(#{call} 1> info_#{name}) 2> err_#{name}"
      end
      puts "Run #{call}"
      system(call)
    else
      puts "Compile first version #{version}"
    end
    return name
  rescue
    print_usage
  end
end

def clean_infos(cleanall = false)
  dirtyfiles = Dir.entries(Dir.pwd).select{|f| f=~/^RAxML_/ or f=~/^err_/ or f=~/^info_/}
  dirtyfiles.delete_if{|f| f =~ /^RAxML_[results|info|states]/} unless cleanall
  dirtyfiles.each{|f| inform_and_delete_file(f)}
end

# tasks -----------------------------------------
desc "print usage"
task :usage do
  print_usage
end
task :clean_compile do 
  inform_and_delete_extension_files(".o")
end

task :clean_infos do 
  clean_infos
end
task :clean_all_infos do 
  clean_infos(cleanall = true)
end

desc "clean all but results"
task :clean => [:clean_compile, :clean_infos]
desc "clean all"
task :cleanall => [:clean_compile, :clean_all_infos]

desc "compile[bayes|std]"
task :compile, [:version] => [:clean_compile] do |t, args|
  args.with_defaults(:version => "bayes")
  inform_and_delete_file(VERSIONS[args.version][:bin]) 
  puts "Compiling version #{args.version}"
  system("make -f #{VERSIONS[args.version][:makefile]}")
end

desc "run[set,version,generations]"
task :run, [:set, :version, :model, :starting_tree, :generations] => [:clean_infos] do |t, args|
  args.with_defaults(:set => "20", :version => "std", :model => "GAMMA", 
                     :generations => 1000, :starting_tree => "treefile")
  puts "run raxml"
  run_name = run_raxml(args.version, args.set, args.model, args.starting_tree, args.generations)
  #puts "  ========  tail info  ============="
  #system("tail -n 22 RAxML_info.#{run_name}")
end

desc "testrun[set,generations]"
task :testrun, [:set, :generations] => [:clean_infos] do |t, args|
  args.with_defaults(:set => "20", :generations => 10000) 
  VERSIONS.keys.each do |version|
    VALID_MODELS.each do |model|
      STARTING_TREE_OPTS.each do |starting_tree|
        run_name = run_raxml(version, args.set, model, starting_tree, args.generations)
      end
    end
  end
end

desc "update ctags"
task :tag do
  system("ctags *.c *.h")
end

task :default => [:usage]
