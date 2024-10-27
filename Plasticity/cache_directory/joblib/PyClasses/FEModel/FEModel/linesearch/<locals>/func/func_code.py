# first line: 1488
        @memory.cache
        def func(alpha,FUNJAC,u,free_ind,h_new,checkSecant=False):
            ux = u.copy()
            ux[free_ind] = u[free_ind] + alpha*h_new
            if checkSecant:
                fb,fc,ft,_,_ = FUNJAC(ux,split=True)
                fb = fb[free_ind]
                fc = fc[free_ind]
                ft = ft[free_ind]
                fb = h_new@fb
                fc = h_new@fc
                ft = h_new@ft

                return fb,fc,ft

            f , _ = FUNJAC(ux)
            f_3 = f[free_ind]
            f3 = h_new@f_3

            return f3, f_3, ux
