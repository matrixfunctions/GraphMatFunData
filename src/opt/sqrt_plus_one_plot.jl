using GraphMatFun, DelimitedFiles
graphnames=["denman_beavers",
            "denman_beavers_opt",
            "ps",
            "ps_opt"]


f=s->sqrt(1+s);
xv=big.(range(-0.5,stop=0.5,length=20));
yv=big.(range(-0.5,stop=0.5,length=20));
for (k,graphname) in enumerate(graphnames)
    graph=import_compgraph("../../data/sqrt_plus_one/$graphname.cgr");

    Z=zeros(size(yv,1),size(xv,1));

    io = open("/home/jarl/jobb/doc/matfun/manuscript/gfx/sqrt_plus_one_surf$k.csv","w");
    for j=1:size(xv,1)
        data=zeros(size(yv,1),3);
        global counter;
        counter=1;
        for i=1:size(yv,1);
            x=xv[i]+1im*yv[j];
            Z[i,j]=abs(eval_graph(graph,x)-f(x));
            Z[i,j]=max(Z[i,j],eps()^2);
            data[counter,:]=[xv[i],yv[j],Z[i,j]];
            counter +=1;
        end
        writedlm(io,data," ");
        write(io,"\n");
    end
    close(io);
end

    #Z[Z .== 0 ] .= eps()/1000;
#plot(z=Z, x=xv, y=yv,Geom.contour(levels=10 .^-float.(1:3:15)))
#plot(z=Z, x=xv, y=yv,Geom.contour(levels=[1e-1;1e-4;1e-7;1e-10;1e-13]))
#writedlm("/tmp/bla.csv",data,' ');
