#!/Users/scott/jbin/Julia-1.0.app/Contents/Resources/julia/bin/julia

###########################################################################
# This program plots the eigenvector loadings from the *.evec file output #
# from Eigenstrat (tested with EIG6.1.4). Any pair of principal components#
# may be plotted provided that the are in the *.evec file.                #
# The input may also include a user file with the sample ids and          #
# populations (e.g. Sample1 Finnish), one sample per line to allow        #
# grouping by population.  A special population may be highlighted in the #
# plot using a fourth argument.                                           #
#                                                                         #
# Requires: Julia 1.0+ , DataFrames, CSV, PyPlot, Colors  and             #
# Anaconda Python is highly recommended for easy PyPlot access.           #
#                                                                         #
# Written by Scott Watkins (2017) GPL                                     #
# Provided to the public in conjunction with Bioinformatics 2018          #
# Revised to julia1.0.3 (2019)                                            #
###########################################################################
println("Hello, thanks for trying plotEIG.jl\n\nLoading packages ...\n\n")


#-------------------- Section: Check user inputs ---------------------#
a = length(ARGS)

if a < 1

    println("USAGE: plotEIG.jl <eigfile.evec> <PCdim1,PCdim2>\n")
    println("Example: plotEIG.jl pca_data.pca.evec 1,2 <pop_codes.txt(opt)> <pop_name(opt)>\n\n    population_codes.txt: a two column file with\n    each sample id and its population name listed on a line.\n\n    pop_name: if used, these samples will be larger in the plot\n    (requires population_codes.txt)\n")
    exit()

end

if a < 2

    println("You must supply the PCs to plot. (e.g. 1,2)\n\nExample: plotEIG.jl pca_data.pca.evec 1,2\n")
    exit();

else

    v = split(ARGS[2], r",")
    println("Plotting data from ", ARGS[1], ": PC ", v[1], " vs PC ", v[2], "\n" )

end

#-------------------- end Check user inputs ---------------------------#

using CSV, DataFrames, Random, PyPlot, Colors


info  = Dict{Any, Any }()         #sample info
pops  = Dict{Any, Any}()          #populations and counts
plist = Any[]                     #pop list
ginfo = Dict{Any,Any}()    
hpops = Dict{Any,Any}()           #special populations


xin = parse(Int64, v[1]) + 1      # cols of df to plot, col1 is names
yin = parse(Int64, v[2]) + 1 

if length(v) == 3
    zin = parse(Int64, v[3]) + 1
end 

eig = CSV.read(ARGS[1], datarow=2, header=false, delim=' ', ignorerepeated=true)
deletecols!(eig, :Column1) 

f=open(ARGS[1])                   #get the data of 2 PCs to plot
a=readlines(f)
close(f)
p = split(a[1],r"\s+")
null = popfirst!(p)


# #------------ Section: assign a color to each population ------------#
if isassigned(ARGS, 3) 
  
    println("Reading population data from ", ARGS[3], "\n")
    open(ARGS[3]) do f    
        for i in eachline(f)
            a = split(i, r"\t|\s+|\,")
            if length(a) > 2
                n = ARGS[3]
                println("File $n should be 2 columns: Sample_id and population_code.\nPlease check the file.\n")
                exit()
            end
            info[a[1]] = a[2]
            ginfo[a[1]] = a[2]
            push!(plist, a[2])	
        end
    end
end

if isassigned(ARGS, 4) 
  
    println("Reading special populations from ", ARGS[4], "\n")
    open(ARGS[4]) do f    
        for i in eachline(f)
            a = split(i, r"\t|\s+|\,")
            if length(a) > 1
                n = ARGS[4]
                println("Filed $n should be a single column with the id of the samples to highlight.\nPlease check the file.\n")
                exit()
            end
            hpops[a[1]] = a[1]
        end
    end
end


#Name the colors for special handling
#green3=[0.0,0.804,0.0]; magenta1=[1.0,0.0,1.0]; dodgerblue3=[0.094,0.455,0.804];
#deepskyblue=[0.1,0.800,1.0]; red2=[0.933,0.0,0.0]; red4=[0.545,0.0,0.0];
#peru=[0.804,0.521,0.247];turquoise=[0.0,0.999,0.999];
#coral1=[1.0,0.447,0.337]; darkorchid1=[0.749,0.243,1.0]; mediumpurple2=[0.537,0.408,0.804];
#goldenrod2=[0.933,0.678,0.055]; greenyellow=[0.678,1.0,0.184];aquamarine1=[0.44,0.999,0.790];
#Random.seed!(9708961)
#colors = shuffle([green3, greenyellow, dodgerblue3, deepskyblue, red2, red4, peru, turquoise, coral1, darkorchid1, mediumpurple2, goldenrod2, magenta1, aquamarine1])

#...OR... Generate colors automatically.
Random.seed!(9703)         #comment to randomize color order
#Seed with white to distiquish against white background, then remove white.
cols = shuffle(distinguishable_colors( length(unique(plist)) + 1, [RGB(1,1,1)] )[2:end])

#cols = distinguishable_colors( length(unique(plist)) )
colors = map(col -> (red(col), green(col), blue(col)), cols)

for c in plist
    pops[c] = get(pops, c, 0) + 1
end

for k in keys(pops)
    println("Population: ", k, " ", pops[k], " samples listed")
    i = popfirst!(colors)
    pops[k] = i              # Each pop gets color code
end

for k in keys(info)          # Transfer pop color to 
    pid = get(info, k, 0)    # each sample to create
    cid = pops[pid]          # color lookup by sample
    info[k] = cid 
end

pop    = Array{String}(undef, size(eig,1))
colors = Array{Any}(undef, size(eig,1)) 

insertcols!(eig, size(eig,2) + 1, :pop => pop)    # Add 2 new cols to eig
insertcols!(eig, size(eig,2) + 1, :color => colors)    # Add 2 new cols to eig
rename!(eig, Dict(:Column2 => :id))

for i in 1:size(eig,1)
    n  = eig[i,:id]              # Add the population and  
    pv = get(ginfo,n,"missing")  # color code corresponding
    pc = get(info, n,"missing")  # to each samples       
    eig[i,:pop] = pv   
    eig[i,:color] = pc	
end

#println(eig)

#------------ End Section: assign a color to each population --------#

if isassigned(ARGS, 3)

   for k in sort(collect(keys(pops)))
       ptable1 = eig[eig[:pop] .== k, :]
       cc = ptable1[:color]
       lab = ptable1[1,:pop]
#      ptable2 = ptable1[ptable1[3] .< 0.043, :]       #user screening if desired 
#      println(ptable2)
       global pointsizes
           if haskey(hpops, ptable1[1,:pop])
               pointsizes=[90]                        #set hpop  to big points
               marker="^"
               edge=[0,0,0]
           else
               pointsizes=[50]                         #set point size for pop plot
               marker="o"
               edge=[0,0,0]
           end  
        plt[:scatter](ptable1[xin], ptable1[yin], s=pointsizes, c=cc, label=lab, alpha=0.6, marker=marker,  edgecolors=edge)
        plt[:legend](loc="center left", bbox_to_anchor=(1, 0.5), markerscale=1.5, labelspacing=1)
    end

else

    n = eig[1]                           #names [:id]            

    x = eig[xin]                         #x values 
    y = eig[yin]                         #y values

    if length(v) == 3
        z = eig[zin]                         #z values if present
    end

    if length(v) == 2
        plt[:scatter](x, y, marker="+", s=30, c="darkslateblue", alpha=0.7, edgecolors="black") 
    end

    if length(v) == 3
#        plt[:scatter](x, y, z,  marker="+", c="darkslateblue", alpha=0.7, edgecolors="black") 
    println("Sorry, 3D plot not yet implemented\n")
    exit
    end

end

plt[:title]("Principal Component Analysis", fontweight="bold", fontsize=16)
plt[:grid](false) 

pcx = v[1]
pcy = v[2]
perx = round(parse(Float64,p[xin]), digits=2)
pery = round(parse(Float64,p[yin]), digits=2)

plt[:xlabel]("PC$pcx\n($perx% variance)", fontsize=14)
plt[:ylabel]("PC$pcy\n($pery% variance)", fontsize=14)

plt[:tight_layout]()               # make plot fit


plt[:show]()                       # plt.show for interactive plot OR

#plt[:savefig]("PCAplot.svg")       # save as vector graphics for editing

#println("PCA plot written to file PCAplot.svg ")
