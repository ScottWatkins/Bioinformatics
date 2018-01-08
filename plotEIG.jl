#!/Users/scott/jbin/Julia-0.6.app/Contents/Resources/julia/bin/julia

###########################################################################
# This program plots the eigenvector loadings from the *.evec file output #
# from Eigenstrat (tested with EIG6.1.4). Any pair of principal components#
# may be plotted provided that the are in the *.evec file.                #
# The input may also include a user file with the sample ids and          #
# populations (e.g. Sample1 Finnish), one sample per line to allow        #
# grouping by population.  A special population may be highlighted in the #
# plot using a fourth argument.                                           #
#                                                                         #
# Requires: Julia 0.6 (+ DataFrames, DataArrays, PyPlot) and              #
# Anaconda Python is highly recommended for easy PyPlot access.           #
#                                                                         #
# Written by Scott Watkins (2017) GPL                                     #
# Provided to the public in conjunction with Bioinformatics 2018          #
###########################################################################

using DataFrames, DataArrays

println("Hello, thanks for trying plotEIG.jl\n")

#-------------------- Section: Check user inputs ---------------------#
if isempty(ARGS)
  println("USAGE: plotEIG.jl <eigfile.evec> <PCdim1,PCdim2>\n")
  println("Example: plotEIG.jl pca_data.pca.evec 1,2 <pop_codes.txt(opt)> <pop_name(opt)>\n\npopulation_codes.txt: a two column file with\neach sample id and its population name listed on a line.\npop_name: if used these samples will be larger in plot (requires population_codes.txt\)\n")
  exit()
end

if isassigned(ARGS, 2)
  v = split(ARGS[2], r",")
  println("Plotting data from ", ARGS[1], ": PC ", v[1], " vs PC ", v[2], "\n" )
else
  println("You must supply the PCs to plot. (e.g. 1,2)")
  exit();
end

if isassigned(ARGS, 4)       # highlight this population 
  Hpop = chomp(ARGS[4])
else
  Hpop = "NULL"
end

#-------------------- end Check user inputs ---------------------------#

info=Dict{Any, Any }()   #sample info
pops = Dict{Any, Any}()  #populations and counts
plist = Any[]            #pop list
ginfo=Dict{Any,Any}()    

xin = parse(Int64, v[1]) + 1      # change input PCs to Int
yin = parse(Int64, v[2]) + 1      # table offset by 1

eig = readtable(ARGS[1], skipstart=1, header=false, separator=' ') #eigenstrat loadings in

f=open(ARGS[1])           #get the variances of 2 PCs to plot
a=readlines(f)
close(f)
p = split(a[1],r"\s+")
null = shift!(p)


#------------ Section: assign a color to each population ------------#
if isassigned(ARGS, 3) 
  
  println("Reading population data from ", ARGS[3], "\n")
  open(ARGS[3]) do f
    for i in eachline(f)
      a = split(i, r"\t|\s+")
      info[a[1]] = a[2]
      ginfo[a[1]] = a[2]
      push!(plist, a[2])	
    end
  end
  #name the colors for easy handling
  green3=[0.0,0.804,0.0]; magenta1=[1.0,0.0,1.0]; dodgerblue3=[0.094,0.455,0.804];
  deepskyblue=[0.1,0.800,1.0]; red2=[0.933,0.0,0.0]; red4=[0.545,0.0,0.0];
  peru=[0.804,0.521,0.247];turquoise=[0.0,0.999,0.999];
  coral1=[1.0,0.447,0.337]; darkorchid1=[0.749,0.243,1.0]; mediumpurple2=[0.537,0.408,0.804];
  goldenrod2=[0.933,0.678,0.055]; greenyellow=[0.678,1.0,0.184];aquamarine1=[0.44,0.999,0.790];
  srand(9708961)
  colors = shuffle([green3, greenyellow, dodgerblue3, deepskyblue, red2, red4, peru, turquoise, coral1, darkorchid1, mediumpurple2, goldenrod2, magenta1, aquamarine1])

  for c in plist
    pops[c]=get(pops,c,0) + 1
  end

  for k in keys(pops)
    println("Population: ", k, " ", pops[k], " samples listed")
    i = shift!(colors)
    pops[k] = i              # Each pop gets color code
  end

  for k in keys(info)        # Transfer pop color to 
    pid = get(info, k, 0)    # each sample to create
    cid = pops[pid]          # color lookup by sample
    info[k] = cid 
  end

  pop = Array{String,1}(size(eig,1))
  colors = Array{Any}(size(eig,1)) 
 
  eig2=hcat(eig, pop, colors)                         # Add 2 new cols to eig
  rename!(eig2, [:x1,:x1_1,:x1_2], [:id,:pop,:color]) # rename id, pop, and color

  r=0                       # populate the table
  for i in eig2[:,:id]      # with the population and  
    pv = get(ginfo,i,"NA")  # color code corresponding
    pc = get(info, i,"NA")  # to each samples
    r=r+1                   # ARGS[3] may contain non-match values
    eig2[r,:pop]=pv   
    eig2[r,:color]=pc	
  end

end
#------------ End Section: assign a color to each population --------#



import PyPlot
const plt = PyPlot

if isassigned(ARGS, 3)

  for k in sort(collect(keys(pops)))
    ptable1 = eig2[eig2[:pop] .== k, :]
    cc = ptable1[:color]
    lab = ptable1[1,:pop]
#    ptable2 = ptable1[ptable1[3] .< 0.043, :]       #user screening if desired 
#    println(ptable2)
    if ptable1[1,:pop] .== Hpop
      pointsizes=[60]                                #set Hpop ARGS[4] to big points
      marker="^"
      edge=[0,0,0]
    else
      pointsizes=[30]                                #set point size for pop plot
      marker="o"
      edge=[0,0,0]
    end  

    plt.scatter(ptable1[xin], ptable1[yin], s=pointsizes, c=cc, label=lab, alpha=0.99, marker=marker,  edgecolors=edge)
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), markerscale=1.5)
  end

else

  n = eig[1]                           #names [:id]            
  x = eig[xin]                         #x values 
  y = eig[yin]                         #y values

  ptsize=[14]
  pointsizes=repmat(ptsize, length(x)) #make array for each pt
  
  for i in 1:length(x)                 # plot the points, ps is point size
    ps=shift!(pointsizes)              # cc is the rbg color code
    cc = get(info, n[i], "Error: sample not found.")
    plt.scatter(x[i], y[i], s=ps, c=[0,0,0]) 
  end

end

plt.title("Principal Component Analysis")
plt.grid("off")                        #grid

pcx = v[1]
pcy = v[2]
perx = round(parse(Float64,p[xin]),1)
pery = round(parse(Float64,p[yin]),1)

plt.xlabel("PC$pcx\n($perx% variance)")
plt.ylabel("PC$pcy\n($pery% variance)")

plt.tight_layout()              # make plot fit
#plt.show()                     # plt.show for interactive plot OR
plt.savefig("plot.svg")         # save as vector graphics for editing
println("PCA plot written to plot.svg.")
