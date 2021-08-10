function degopt=read_degopt(fname,use_vpa)
% Usage
%   degopt=read_degopt(fname,[true|false])

%fname = '../data/exp/exp_sastre_m5_opt_rho1_9.cgr';
%fname = '../data/exp/exp_ps_m8_opt_rho13_5.cgr';

    fid = fopen(fname);
    rows={};
    n=50;
    maxrow=0;
    Ha=zeros(n,n);
    Hb=zeros(n,n);
    y=zeros(n+2,1);

    %    use_vpa=true;
    if (use_vpa)
        % Use high precision
        convertnum=@(s) vpa(s);
        Ha=vpa(Ha);
        Hb=vpa(Hb);
        y=vpa(y);
    else
        convertnum=@(s) str2num(s);
    end

    while true
        thisline = fgetl(fid);
        if ~ischar(thisline); break; end  %end of file
                                          %now check whether the string in thisline is a "word", and store it if it is.
                                          %then

        x=regexp(thisline,'[=;]','split');
        if (length(x)>1)
            lhs=x{1};
            rhs=x{2};
            if (startsWith(thisline,'coeff1'))
                coeff1=convertnum(rhs);
            end
            if (startsWith(thisline,'coeff2'))
                coeff2=convertnum(rhs);
            end

            if (startsWith(rhs,'coeff1'));
                % It's a linear combination
                idxs=[];
                xx=regexp(lhs,'B[ab](\d+)_(\d+)','tokens');
                if (length(xx)>0)
                    xx1=xx{1};
                    row=str2num(xx1{1})-1;
                    z=str2num(xx1{2});
                    if (z == 2)
                        idxs=[1 2];
                        vals=[coeff1 coeff2];
                    else
                        idxs=[z];
                        vals=[coeff2];
                    end
                else
                    xx=regexp(lhs,'B[ab](\d+)','tokens');
                    if (length(xx)>0)
                        xx1=xx{1};
                        z=str2num(xx1{1});

                        row=z-1;
                        if (z==2)
                            idxs=[1 2];
                            vals=[coeff1 coeff2];
                        else
                            idxs=[z];
                            vals=[coeff2];
                        end

                    end
                end

                for i=1:length(idxs)
                    j=idxs(i);
                    v=vals(i);
                    if (startsWith(lhs,'Ba'))
                        Ha(row,j)=v;
                    else
                        Hb(row,j)=v;
                    end
                end

                maxrow=max(row,maxrow);

                %            if (z==2)
                %                idxs
                %                vals
                %                asd
                %            end


            end


            if (startsWith(lhs,'T2k'));
                xx=regexp(lhs,'T2k(\d+)','tokens');
                xx1=xx{1};
                z=str2num(xx1{1});
                if (z==2)
                    y(1)=coeff1;
                    y(2)=coeff2;
                else
                    y(z)=coeff2;
                end
            end


        end
    end
    fclose(fid);

    Ha=Ha(1:maxrow,1:maxrow+1)
    Hb=Hb(1:maxrow,1:maxrow+1)
    y(maxrow+2)=y(maxrow+3);
    y=y(1:maxrow+2);

    degopt=struct('Ha',Ha,'Hb',Hb,'y',y);
end
