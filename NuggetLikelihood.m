function[val] = NuggetLikelihood(x,Phi_Phi,Phi_S_Phi,trS,n)

        l = size(Phi_S_Phi,1);
        
        Phi_Phi_over_nug = Phi_Phi/x(1);
        Phi_S_Phi_over_nugsq = Phi_S_Phi/(x(1)^2);

        for j=1:l
            Phi_Phi_over_nug(:,j) =   Phi_Phi_over_nug(:,j)*1/sqrt(x(2)+1/(x(1)));
            Phi_S_Phi_over_nugsq(:,j) = Phi_S_Phi_over_nugsq(:,j) * 1/(x(2)+1/(x(1)));
        end
        cholQ_plus_Phi = chol(x(2)*eye(l)+Phi_Phi_over_nug);
        %Working on making this work with nonorthogonal basis fns
        logdetpart = 2*sum(log(diag(cholQ_plus_Phi))) - l*log(x(2));
        tracepart  = norm(cholQ_plus_Phi'\Phi_S_Phi_over_nugsq,'fro')^2/l;

        constantpart= n*log(x(1)) + trS/(x(1));

        
        val =  logdetpart - tracepart + constantpart;
end

