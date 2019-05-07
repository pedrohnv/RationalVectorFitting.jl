using LinearAlgebra

function flag_poles(poles)
    """ Finding out which starting poles are complex """
    N = length(poles);
    cindex = zeros(Int, N);
    for m=1:N
        if abs(imag(poles[m])) >= eps(Float64)
            if m==1
                cindex[m] = 1;
            else
                if cindex[m-1]==0 || cindex[m-1]==2
                    cindex[m] = 1;
                    cindex[m+1] = 2;
                else
                    cindex[m] = 2;
                end
            end
        end
    end
    return cindex
end

function vectfit3(f, s, poles, weight, relax=true, stable=true, asymp=2,
                  skip_pole=false, skip_res=false)
    """
         ===========================================================
         =   Fast Relaxed Vector Fitting                           =
         =   Version 1.0                                           =
         =   Last revised: 08.08.2008                              =
         =   Written by: Bjorn Gustavsen                           =
         =   SINTEF Energy Research, N-7465 Trondheim, NORWAY      =
         =   bjorn.gustavsen@sintef.no                             =
         =   http://www.energy.sintef.no/Produkt/VECTFIT/index.asp =
         =   Note: RESTRICTED to NON-COMMERCIAL use                =
         ===========================================================

     PURPOSE : Approximate f(s) with a state-space model

             f(s)=C*(s*I-A)^(-1)*B +D +s*E

               where f(s) is a singe element or a vector of elements.
               When f(s) is a vector, all elements become fitted with a common
               pole set.

     INPUT :

     f(s) : function (vector) to be fitted.
            dimension : (Nc,Ns)
                         Nc : number of elements in vector
                         Ns : number of frequency samples

     s : vector of frequency points [rad/sec]
            dimension : (1,Ns)

     poles : vector of initial poles [rad/sec]
             dimension : (1,N)

     weight: the rows in the system matrix are weighted using this array. Can be used
             for achieving higher accuracy at desired frequency samples.
             If no weighting is desired, use unitary weights: weight=ones(1,Ns).

             Two dimensions are allowed:
               dimension : (1,Ns) --> Common weighting for all vector elements.
               dimension : (Nc,Ns)--> Individual weighting for vector elements.

     relax: optional, default=true
        true --> Use relaxed nontriviality constraint
        false --> Use nontriviality constraint of "standard" vector fitting

     stable: optional, default=true
         false --> unstable poles are kept unchanged
         true --> unstable poles are made stable by 'flipping' them
                       into the left half-plane

     asymp: optional, default=2
         1 --> Fitting with D=0,  E=0
         2 --> Fitting with D!=0, E=0
         3 --> Fitting with D!=0, E!=0

     skip_pole: optional, default=false
         true --> The pole identification part is skipped, i.e (C,D,E)
                   are identified using the initial poles (A) as final poles.

     skip_res: optional, default=false
         true --> The residue identification part is skipped, i.e. only the
                    poles (A) are identified while C,D,E are returned as zero.
     ========
     OUTPUT :

         fit(s) = C*(s*I-(A)^(-1)*B +D +s.*E

     poles(1,N)    : new poles

     rmserr(1) : root-mean-square error of approximation for f(s).
                       (0 is returned if skip_res==1)
     fit(Nc,Ns): Rational approximation at samples. (0 is returned if
                 skip_res==1).

     APPROACH:
     The identification is done using the pole relocating method known as Vector Fitting [1],
     with relaxed non-triviality constraint for faster convergence and smaller fitting errors [2],
     and utilization of matrix structure for fast solution of the pole identifion step [3].

    ********************************************************************************
     NOTE: The use of this program is limited to NON-COMMERCIAL usage only.
     If the program code (or a modified version) is used in a scientific work,
     then reference should be made to the following:

     [1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency
         domain responses by Vector Fitting", IEEE Trans. Power Delivery,
         vol. 14, no. 3, pp. 1052-1061, July 1999.

     [2] B. Gustavsen, "Improving the pole relocating properties of vector
         fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592,
         July 2006.

     [3] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter,
         "Macromodeling of Multiport Systems Using a Fast Implementation of
         the Vector Fitting Method", IEEE Microwave and Wireless Components
         Letters, vol. 18, no. 6, pp. 383-385, June 2008.
    ********************************************************************************
     This example script is part of the vector fitting package (VFIT3.zip)
     Last revised: 08.08.2008.
     Created by:   Bjorn Gustavsen.
    """
    #Tolerances used by relaxed version of vector fitting
    TOLlow = 1e-18;
    TOLhigh = 1e18;
    a = 0;
    b = 0;
    try
        a, b = size(poles);
    catch
        a = 1;
        b = size(poles)[1];
    end
    if s[1]==0 && a==1
        if poles[1]==0 && poles[2]!=0
            poles[1] = -1;
        elseif poles[2]==0 && poles[1]!=0
            poles[2] = -1;
        elseif poles[1]==0 && poles[2]==0
            poles[1] = -1+10im;
            poles[2] = -1-10im;
        end
    end
    rmserr = [];
    #try
    #    a, b = size(s);
    #catch
    #    a = 1;
    #    b = size(s)[1];
    #end
    #if a < b
        s = transpose(s);
    #end
    #TODO sanity check the arguments
    LAMBD = diagm(0 => poles);
    Ns = length(s);
    N = size(LAMBD)[1]; #length(LAMBD);
    Nc = length(f[:,1]);
    B = ones(N,1);
    SERA = poles;
    SERC = zeros(Nc,N);
    SERD = zeros(Nc,1);
    SERE = zeros(Nc,1);
    roetter = poles;
    fit = zeros(ComplexF64, Nc, Ns);

    weight = transpose(weight);
    if length(weight[1,:])==1
        common_weight = 1;
    elseif length(weight[1,:])==Nc
        common_weight = 0;
    else
        display("ERROR in vectfit3.m: Invalid size of array weight")
        return
    end

    if asymp==1
        offs = 0;
    elseif asymp==2
        offs = 1;
    else
        offs = 2;
    end
    # =======================================================
    ## POLE IDENTIFICATION:
    # =======================================================
    if !(skip_pole)
        Escale = zeros(1,Nc+1);
        # =======================================================
        # Finding out which starting poles are complex :
        # =======================================================
        cindex = flag_poles(poles);
        # =======================================================
        ## Building system - matrix :
        # =======================================================
        if asymp==1 || asymp==2
            Dk = zeros(ComplexF64, Ns,N+1);
            for i=1:Ns
                Dk[i,N+1] = 1;
            end
        elseif asymp==3
            Dk = zeros(ComplexF64, Ns,N+2);
            for i=1:Ns
                Dk[i,N+1] = 1;
                Dk[i,N+2] = s[i];
            end
        else
            Dk = zeros(ComplexF64, Ns,N);
        end
        for m=1:N
            if cindex[m]==0      #real pole
                Dk[:,m] = 1 ./(s .- LAMBD[m,m]);
            elseif cindex[m]==1  #complex pole, 1st part
                Dk[:,m] = 1 ./(s .- LAMBD[m,m]) + 1 ./(s .- conj(LAMBD[m,m]));
                Dk[:,m+1] = 1im ./(s .- LAMBD[m,m]) - 1im ./(s .- conj(LAMBD[m,m]));
            end
        end

        #Scaling for last row of LS-problem (pole identification)
        scale = 0;
        for m=1:Nc
            if length(weight[1,:])==1
                scale = scale + (norm(weight[:,1].*transpose(f[m,:])))^2;
            else
                scale = scale + (norm(weight[:,m].*transpose(f[m,:])))^2;
            end
        end
        scale = sqrt(scale)/Ns;

        if relax
            #Escale=zeros(1,Nc*(N+offs)+N+1);
            #scale=norm(f);#/Ns; #Scaling for sigma in LS problem
            AA = zeros(ComplexF64, Nc*(N+1), N+1);
            bb = zeros(ComplexF64, Nc*(N+1), 1);
            Escale = zeros(ComplexF64, 1, length(AA[1,:]));

            for n=1:Nc
                A = zeros(ComplexF64, Ns, (N+offs)+N+1); #b=zeros(Ns*Nc+1,1);
                if common_weight==1
                    weig = weight;
                else
                    weig = weight[:,n];
                end
                for m=1:N+offs #left block
                    A[1:Ns,m] = weig .* Dk[1:Ns,m];
                end
                inda = N + offs;
                for m=1:N+1 #right block
                    A[1:Ns,inda+m] = -weig .* Dk[1:Ns,m] .* f[n,1:Ns];
                end
                A = [real(A); imag(A); zeros(1, (N+offs)+N+1)];
                #Integral criterion for sigma:
                offset = (N+offs);
                if n==Nc
                    for mm=1:N+1
                        A[2*Ns+1, offset+mm] = real(scale*sum(Dk[:,mm]));
                    end
                end
                F = qr(A);
                Q = Matrix(F.Q);
                R = Matrix(F.R);
                ind1 = N + offs + 1;
                ind2 = N + offs + N + 1;
                R22 = R[ind1:ind2,ind1:ind2];
                AA[(n-1)*(N+1)+1:n*(N+1),:] = R22;
                if n==Nc
                    bb[(n-1)*(N+1)+1:n*(N+1),1] = (Q[end,N+offs+1:end]')*Ns*scale;
                end
            end #for n=1:Nc
            for col=1:length(AA[1,:])
                Escale[col] = 1/norm(AA[:,col]);
                AA[:,col] = Escale[col] .* AA[:,col];
            end
            x = AA\bb;
            #size(x),size(Escale)
            x = x .* transpose(Escale);
        end #if relax

        #Situation: No relaxation, or produced D of sigma extremely small and large. Solve again, without relaxation
        if !relax | (abs(x[end]) < TOLlow) | (abs(x[end]) > TOLhigh)
            AA = zeros(ComplexF64, Nc*(N),N);
            bb = zeros(ComplexF64, Nc*(N),1);
            if !relax
                Dnew = 1;
            else
                if x[end]==0
                    Dnew = 1;
                elseif abs(x[end])<TOLlow
                    Dnew = sign(x[end])*TOLlow;
                elseif abs(x[end])>TOLhigh
                    Dnew = sign(x[end])*TOLhigh;
                end
            end

            for n=1:Nc
                A=zeros(Ns,(N+offs) +N); #b=zeros(Ns*Nc+1,1);
                Escale=zeros(1,N);
                if common_weight==1
                    weig = weight;
                else
                    weig=weight[:,n];
                end

                for m=1:N+offs #left block
                    A[1:Ns,m] = weig.*Dk[1:Ns,m];
                end
                inda = N+offs;
                for m=1:N #right block
                    A[1:Ns,inda+m] = -weig.*Dk[1:Ns,m].*transpose(f[n,1:Ns]);
                end
                b = Dnew*weig.*transpose(f[n,1:Ns]);
                A = [real(A);imag(A)];
                b = [real(b);imag(b)];
                offset = (N+offs);
                Q, R = qr(A,0);
                ind1 = N+offs+1;
                ind2 = N+offs+N;
                R22 = R[ind1:ind2,ind1:ind2];
                AA[(n-1)*N+1:n*N,:] = R22;
                bb[((n-1)*N+1:n*N,1)] = transpose(Q[:,ind1:ind2])*b;
            end #for n=1:Nc
            for col=1:length(AA[1,:])
                Escale[col] = 1 ./ norm(AA[:,col]);
                AA[:,col] = Escale[col] .* AA[:,col];
            end
            #if use_normal==1
            #  x=AA.'*AA\(AA.'*bb);
            #else
            x = AA\bb;
            #end
            x = x .* transpose(Escale);
            x = [x;Dnew];
        end  #if relax==0 | abs(x[end])<TOLlow | abs(x[end])>TOLhigh
        # ************************************
        C = x[1:end-1];
        D = x[end]; #NEW!!
        #We now change back to make C complex :
        # **************
        for m=1:N
            if cindex[m]==1
                for n=1:1 #Nc+1
                    r1 = C[m];
                    r2 = C[m+1];
                    C[m] = r1 + 1im*r2;
                    C[m+1] = r1 - 1im*r2;
                end
            end
        end
        # **************

        # =============================================================================
        # We now calculate the zeros for sigma :
        # =============================================================================
        #oldLAMBD=LAMBD;oldB=B;oldC=C;7
        m = 0;
        for n=1:N
            m = m+1;
            if m<N
                if ( abs(LAMBD[m,m]) > abs(real(LAMBD[m,m])) ) #complex number?
                    LAMBD[m+1,m] = -imag(LAMBD[m,m]);
                    LAMBD[m,m+1] = imag(LAMBD[m,m]);
                    LAMBD[m,m] = real(LAMBD[m,m]);
                    LAMBD[m+1,m+1] = LAMBD[m,m];
                    B[m,1] = 2;
                    B[m+1,1] = 0;
                    koko = C[m];
                    C[m] = real(koko);
                    C[m+1] = imag(koko);
                    m = m+1;
                end
            end
        end
        ZER = LAMBD - B*transpose(C)/D;
        ZER = Array{Float64,2}(ZER); #FIXME that should not be needed, but gives non-conjugate complex pairs...
        roetter = transpose(eigvals(ZER));
        unstables = map(r -> real(r) > 0, roetter);
        if stable
            roetter[unstables] = roetter[unstables] - 2*real(roetter[unstables]); #Forcing unstable poles to be stable...
        end
        roetter = transpose(roetter);
        sort!(roetter, by = x -> (real(x), imag(x)));
        N = length(roetter);
        # =============================================
        #Sorterer polene s.a. de reelle kommer først:
        for n=1:N
            for m=n+1:N
                if (abs(imag(roetter[m])) <= eps(Float64)
                    && abs(imag(roetter[n])) >= eps(Float64))
                    trans = roetter[n];
                    roetter[n] = roetter[m];
                    roetter[m] = trans;
                end
            end
        end
        N1=0;
        for m=1:N
            if imag(roetter[m])==0
                N1 = m;
            end
        end
        if N1<N
            #roetter[N1+1:N] = sort(roetter[N1+1:N]);
            sort!(roetter[N1+1:N], by = x -> (real(x), imag(x)));
        end       # N1: n.o. real poles
        #N2=N-N1; # N2: n.o. imag.poles
        roetter = roetter - 2im*imag(roetter); #10.11.97 !!!
        SERA = transpose(roetter);
    end #if !skip_pole
    # =========================================================================
    #  RESIDUE IDENTIFICATION:
    # =========================================================================
    if !(skip_res)
        # ============================================================================
        # We now calculate SER for f, using the modified zeros of sigma as new poles :
        # ============================================================================
        LAMBD = diagm(0 => roetter);
        #B=ones(N,1);
        # Finding out which poles are complex :
        cindex = flag_poles(roetter);

        # ===============================================================================
        # We now calculate the SER for f (new fitting), using the above calculated
        # zeros as known poles :
        # ===============================================================================
        if asymp==1
            A = zeros(ComplexF64, 2*Ns,N);
            BB = zeros(ComplexF64, 2*Ns,Nc);
        elseif asymp==2
            A = zeros(ComplexF64, 2*Ns,N+1);
            BB = zeros(ComplexF64, 2*Ns,Nc);
        else
            A = zeros(ComplexF64, 2*Ns,N+2);
            BB = zeros(ComplexF64, 2*Ns,Nc);
        end
        for m=1:N
            if cindex[m]==0      #real pole
                Dk[:,m] = 1 ./(s .- LAMBD[m,m]);
            elseif cindex[m]==1  #complex pole, 1st part
                Dk[:,m] = 1 ./(s .- LAMBD[m,m]) + 1 ./(s .- conj(LAMBD[m,m]));
                Dk[:,m+1] = 1im ./(s .- LAMBD[m,m]) - 1im ./(s .- conj(LAMBD[m,m]));
            end
        end
        if common_weight==1
            weight = transpose(weight);
            #I3=diag(ones(1,Nc));I3(:,Nc)=[];
            Dk = zeros(ComplexF64, Ns,N);
            for m=1:N
                if cindex[m]==0      #real pole
                    Dk[:,m] = weight ./(s .- LAMBD[m,m]);
                elseif cindex[m]==1  #complex pole, 1st part
                    Dk[:,m] = weight ./(s .- LAMBD[m,m]) + weight ./(s .- conj(LAMBD[m,m]));
                    Dk[:,m+1] = 1im .*weight ./(s .- LAMBD[m,m]) - 1im .*weight ./(s .- conj(LAMBD[m,m]));
                end
            end
            if asymp==1
                A[1:Ns,1:N] = Dk;
            elseif asymp==2
                A[1:Ns,1:N] = Dk;
                A[1:Ns,N+1] = weight;
            else
                A[1:Ns,1:N] = Dk;
                A[1:Ns,N+1] = weight;
                A[1:Ns,N+2] = weight.*s;
            end
            for m=1:Nc
                BB[1:Ns,m] = weight.*transpose(f[m,:]);
            end
            A[Ns+1:2*Ns,:] = imag(A[1:Ns,:]);
            A[1:Ns,:] = real(A[1:Ns,:]);
            BB[Ns+1:2*Ns,:] = imag(BB[1:Ns,:]);
            BB[1:Ns,:] = real(BB[1:Ns,:]);

            if asymp==2
                A[1:Ns,N+1] = A[1:Ns,N+1];
            elseif asymp==3
                A[1:Ns,N+1] = A[1:Ns,N+1];
                A[Ns+1:2*Ns,N+2] = A[Ns+1:2*Ns,N+2];
            end
            Escale = zeros(ComplexF64, 1,length(A[1,:]));
            for col=1:length(A[1,:]);
                Escale[col] = norm(A[:,col],2);
                A[:,col] = A[:,col]./Escale[col];
            end
            X = A\BB;
            for n=1:Nc
                X[:,n] = X[:,n]./transpose(Escale);
            end
            X = transpose(X);
            C = X[:,1:N];

            if asymp==2
                SERD = X[:,N+1];
            elseif asymp==3
                SERE = X[:,N+2];
                SERD = X[:,N+1];
            end

        else #if common_weight==1
            SERD = zeros(ComplexF64, Nc,1);
            SERE = zeros(ComplexF64, Nc,1);
            C = zeros(ComplexF64, Nc,N);
            for n=1:Nc
                if asymp==1
                    A[1:Ns,1:N] = Dk;
                elseif asymp==2
                    A[1:Ns,1:N] = Dk;
                    A[1:Ns,N+1] = 1;
                else
                    A[1:Ns,1:N] = Dk;
                    A[1:Ns,N+1] = 1;
                    A[1:Ns,N+2] = s;
                end
                for m=1:length(A[1,:])
                    A[1:Ns,m] = weight[:,n].*A[1:Ns,m];
                end
                BB = weight[:,n].*transpose(f(n,:));
                A[Ns+1:2*Ns,:] = imag(A[1:Ns,:]);
                A[1:Ns,:] = real(A[1:Ns,:]);
                BB[Ns+1:2*Ns] = imag(BB[1:Ns]);
                BB[1:Ns] = real(BB[1:Ns]);

                if asymp==2
                    A[1:Ns,N+1] = A[1:Ns,N+1];
                elseif asymp==3
                    A[1:Ns,N+1] = A[1:Ns,N+1];
                    A[Ns+1:2*Ns,N+2] = A[Ns+1:2*Ns,N+2];
                end
                Escale = zeros(1,length(A[1,:]));
                for col=1:length(A[1,:]);
                    Escale[col] = norm(A[:,col],2);
                    A[:,col] = A[:,col]./Escale[col];
                end
                x = A\BB;
                x = x./transpose(Escale);
                C[n,1:N] = transpose(x[1:N]);
                if asymp==2
                    SERD[n] = x[N+1];
                elseif asymp==3
                    SERE[n] = x[N+2];
                    SERD[n] = x[N+1];
                end
            end #for n=1:Nc
        end #if common_weight==1
        # =========================================================================
        #We now change back to make C complex.
        for m=1:N
            if cindex[m]==1
                for n=1:Nc
                    r1 = C[n,m];
                    r2 = C[n,m+1];
                    C[n,m] = r1 + 1im*r2;
                    C[n,m+1] = r1 - 1im*r2;
                end
            end
        end
        # **************
        B = ones(N,1);
        # ====================================================
        SERA = LAMBD;
        SERB = B;
        SERC = C;
        Dk = zeros(ComplexF64, Ns,N);
        for m=1:N
            Dk[:,m] = 1 ./(s .- SERA[m,m]);
        end
        for n=1:Nc
            fit[n,:] = transpose(Dk*SERC[n,:]);
            if asymp==2
                fit[n,:] = fit[n,:] .+ SERD[n];
            elseif asymp==3
                fit[n,:] = fit[n,:] .+ SERD[n] .+ transpose(s).*SERE[n];
            end
        end
        #fit = transpose(fit);
        #f = transpose(f);
        dif1 = fit - f;
        rmserr = sqrt( sum( sum(map(abs, dif1.^2)) ) )/sqrt(Nc*Ns);
        #fit = transpose(fit);
    end # if !skip_res
    A = SERA;
    poles = A;
    if !(skip_res)
        B = SERB;
        C = SERC;
        D = SERD;
        E = SERE;
    else
        B = ones(N,1);
        C = zeros(Nc,N);
        D = zeros(Nc,Nc);
        E = zeros(Nc,Nc);
        rmserr = 0;
    end;
    #A = sparse(diag(A));    % A is complex, make it diagonal
    #SER.A=A; SER.B=B; SER.C=C; SER.D=D; SER.E=E;
    return diag(poles), rmserr, fit
end
