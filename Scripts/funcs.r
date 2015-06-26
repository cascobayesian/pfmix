run.mcmc = function(num.test,allele.freq,data.set,num.iter,thin)
{
	model.mcmc = rep(list(),num.test)

	cat(paste("\n\nRunning MCMC for 1-",num.test," components on ",sample.name,"...\n\n",sep=""))
	for(j in 1:num.test)
	{
		model.mcmc[[j]] = mh.mcmc(j, allele.freq, data.set, num.iter,thin)
	}

	min.val = 1e10;
	elt = 0;
	bic = rep(0,num.test)

	for (j in 1:num.test)
	{
		max.llk = find.max.llk(model.mcmc[[j]])
		bic[j] = as.numeric(-2*max.llk+(j+2)*(log(dim(data.set)[1])));
		elt = which.min(bic);

	}
	model.out = model.mcmc[[elt]];
}

plot.figure = function(current,sample.name,data.set,allele.freq)
{
	
	figure.file = paste("Figures/sample",sample.name,"pdf",sep=".")
	#print(figure.file)
	pdf(figure.file)
	current = current[[1]][[1]]
	p.mat = generate.p.mat(current$p,current$alpha,allele.freq);
	l     = current$lambda;
	num.comp = current$num.comp

	num.elt = 2^num.comp

	par(cex.axis=1.3,cex.lab=1.5,mar=c(5,5,4,2)+0.1)
	for (i in 1:num.elt)
	{
		if (i == 1)
		{
			plot(allele.freq,p.mat[,i],xlim=c(-0.01,1.01),ylim=c(-0.01,1.01),col=grey(1-l[,i]),pch=20,lwd=0.5,xlab="Population-level allele frequency",ylab="Within-sample allele frequency",xaxs="i",yaxs="i");
		}
		else
		{
			points(allele.freq,p.mat[,i],col=grey(1-l[,i]),pch=20,lwd=0.5);
		}
	}

	points(allele.freq,p.mat[,i],xlim=c(0,1),ylim=c(0,1),col=grey(1-l[,i]),pch=20,lwd=0.5);
	points(allele.freq,data.set[,1]/data.set[,2],xlim=c(0,1),ylim=c(0,1),col=rgb(0.1,0.1,1),pch=20,cex=1);

	text(x=0.1,y=0.9,labels=paste("K=",current$num.comp,sep=""),cex=2,col=grey(0.1))
	for(i in 1:current$num.comp)
	{
		text(x=0.15+(i-1)*0.1,y=0.8,labels=round(sort(current$p[2:(2+current$num.comp)]),2)[i],cex=1.5,col=grey(0.2))
	}

	dev.off()
}


find.max.llk = function(chain)
{
	num.iter = length(chain)
	llk = rep(NA,num.iter)
	for(i in 1:num.iter)
	{
		llk[i] = chain[[i]][[2]]
	}
	return(max(llk))
}

generate.lambda = function(num.comp, allele.freq)
{
	num.col = 2^num.comp;
	lambda = matrix(rep(0,num.col*length(allele.freq)),ncol=num.col)

	lambda.comp = matrix(rep(0,(num.comp+1)*length(allele.freq)),ncol=(num.comp+1))

	total = 0*allele.freq;
	for (i in 0:num.comp)
	{
		lambda.comp[,i+1] = allele.freq^i*(1-allele.freq)^(num.comp-i);
		bounds = calc.bounds(num.comp,i);
		lambda[,bounds[1]:bounds[2]] = lambda.comp[,i+1];
	}
	return(lambda)
}

generate.nu = function(num.comp, allele.freq)
{
	num.col = 2^num.comp;
	lambda = matrix(rep(0,num.col*length(allele.freq)),ncol=num.col)

	lambda.comp = matrix(rep(0,(num.comp+1)*length(allele.freq)),ncol=(num.comp+1))

	total = 0*allele.freq;
	for (i in 0:num.comp)
	{
		lambda.comp[,i+1] = allele.freq^i*(1-allele.freq)^(num.comp-i);
		bounds = calc.bounds(num.comp,i);
		lambda[,bounds[1]:bounds[2]] = lambda.comp[,i+1];
	}
	return(lambda)
}

calc.llk = function(data.set,model,allele.freq)
{
	num.snps = dim(data.set)[1]
	num.comp = model$num.comp

	num.ent = 2^num.comp

	val = matrix(rep(0,num.snps*num.ent),ncol=num.ent)

	for(i in 1:num.ent)
	{	val[,i] = model$lambda[,i]*dbeta.binom.zi(data.set[,1],data.set[,2],model$p.mat[,i],model$nu)	}	

	return(sum(log(rowSums(val))));
}

draw = function()
{
	## l : lambda, a : alpha, p : ps, r : paras, i : pi, n
	d = sample(c('l','a','p','r','i','n'),1,prob=c(0,1,1,0,0,1));
	return(d);
}


generate.p = function(num.comp)
{
	p = t(rdiric(1,rep(1,num.comp)));		
	ps = generate.subsets(num.comp);
	p.e = rep(0,0)
	for (i in 1:(2^num.comp))
	{				
		p.e = c(p.e,sum(p[ps[[i]]]))
	}
	return(p.e)
}

generate.p.merge = function(model)
{
	p = model$p;
	num.comp = model$num.comp
	if (num.comp>1)
	{
		these = sample(1:num.comp,2,replace=FALSE)
		total = sum(p[these+1])
		u     = runif(1);
		p[these+1] = c(total*u,total*(1-u));
	}
	ps = generate.subsets(num.comp);
	p.d = p[(2):(num.comp+1)];
	p.e = rep(0,0)
	for (i in 1:(2^num.comp))
	{				
		p.e = c(p.e,sum(p.d[ps[[i]]]))
	}
	return(p.e)
}

generate.p.alpha = function(q,alpha,allele.freq)
{
	return(q*(1-alpha)+alpha*allele.freq);
}

generate.p.mat = function(p,alpha,allele.freq)
{
	return(sapply(p,generate.p.alpha,alpha,allele.freq))
}



init.simple = function(num.comp,allele.freq, alpha=0.1, nu=25,p=NA,val=1,pi=NA)
{
	if (is.na(p[1]))
	{
		p = generate.p(num.comp)
	}
	pi = c(0,1,0)#t(rdirichlet(1,rep(1,3)));	
	paras = c(1,1)


	p.mat = generate.p.mat(p,alpha,allele.freq);
	lambda = generate.lambda(num.comp,allele.freq)

	model = list(num.comp,p,p.mat,lambda,pi,paras,alpha,nu);
	
	class(model) = "model";
	names(model)  = c("num.comp","p","p.mat","lambda","pi","paras","alpha","nu");
	return(model)
}
	
update.model = function(model,allele.freq)
{
	d = draw();
	if (d == 'a')
	{
		model = update.alpha(model);
	}else if (d == 'l')
	{
		model = update.lambda(model,allele.freq);
	}else if (d == 'r')
	{
		model = update.paras(model);
	}else if (d == 'p')
	{
		model = update.p(model,allele.freq);
	}else if (d == 'i')
	{
#		print(d)
		model = update.pi(model);
	}else if (d == 'n')
	{
		model = update.nu(model);
	}else
	{
		print("No draw!");	
		exit();
	}

	return(model);	
}

update.nu = function(model)
{
	model$nu = 20+rexp(1,0.1);
	return(model)
}


update.lambda= function(model,allele.freq)
{
	model$lambda = generate.lambda(num.comp,allele.freq);
	return(model);
}

update.paras = function(model)
{
	model$paras = rexp(2,1);
	return(model);
}

update.pi = function(model)
{
	model$pi = t(rdiric(1,rep(1,3)));
	return(model);
}

update.alpha = function(model)
{
	model$alpha = runif(1);
	model$p.mat = generate.p.mat(model$p,model$alpha,allele.freq);
	return(model);
}

update.p = function(model,allele.freq)
{
	model$p = generate.p.merge(model);
	model$p.mat = generate.p.mat(model$p,model$alpha,allele.freq)
	return(model);	
}

mle = function(num.comp, allele.freq,data.set, num.iter = 20000,thin=400)
{
	current = init.simple(num.comp,allele.freq)
	proposed = current;

	curr.llk = calc.llk(data.set,current,allele.freq);
	prop.llk = curr.llk

#	outs = rep(list(),num.iter/thin)

	for (i in 1:num.iter)
	{
		proposed = current;
		proposed = update.model(proposed,allele.freq);
		prop.llk = calc.llk(data.set,proposed,allele.freq);

		if (prop.llk > curr.llk - 5*max((num.iter-2*i)/num.iter,0))
		{
			curr.llk = prop.llk;
			current  = proposed;
		}
		if (i %% thin ==0)
		{
			print(c(i,curr.llk))
		}
	}
	return(list(current,curr.llk))	
}

generate.subsets = function(num.comp,i)
{
	set = powerset(c(1:num.comp));
	class(set) = "ps";
	set = sort(set);
	return(set)
}

powerset = function(s){
    len = length(s)
    l = vector(mode="list",length=2^len) ; l[[1]]=numeric()
    counter = 1L
    for(x in 1L:length(s)){
        for(subset in 1L:counter){
            counter=counter+1L
            l[[counter]] = c(l[[subset]],s[x])
        }
    }
    return(l)
}

calc.bounds = function(num.comp,i)
{
	if (i == 0)
	{	return(c(1,1));	
	}else if (num.comp == i)
	{
		return(c(2^num.comp,2^num.comp));
	}else
	{
		fore = sum(choose(num.comp,0:(i-1)))+1;
		aft  = sum(choose(num.comp,0:i));
		return(c(fore,aft));
	}
}

### initialize.model

init.model = function(num.comp,allele.freq, alpha=0.1, nu=25,p=NA,val=1)
{
	if (is.na(p[1]))
	{	p = generate.p(num.comp)	}

	p.mat  = generate.p.mat(p,alpha,allele.freq);
	lambda = generate.lambda(num.comp,allele.freq)
	model  = list(num.comp,p,p.mat,lambda,alpha,nu);
	
	class(model)  = "model";
	names(model)  = c("num.comp","p","p.mat","lambda","alpha","nu");
	return(model)
}

### draw function

draw.mcmc = function()
{
	## l : lambda, a : alpha, p : ps, r : paras, i : pi, v : val, n
	d = sample(c('a','p','n'),1,prob=c(2,4,1));
	return(d);
}


### calc.prior

calc.prior = function(model)
{
	log.prior.alpha = 0;
	log.prior.p     = 0;
	log.prior.nu    = dexp(model$nu,0.1,log=TRUE)	
	return(sum(log.prior.alpha + log.prior.p + log.prior.nu));	
}

### calc.prop

calc.proposal = function(model.new,model.old)
{
	log.q.alpha = 0;
	log.q.p     = 0;
	log.q.nu    = dexp(model.old$nu,0.1,log=TRUE) - dexp(model.new$nu,0.1,log=TRUE); 
	return(log.q.alpha + log.q.p + log.q.nu);
}

### update model

update.model.mcmc = function(model)
{
	d = draw.mcmc();
	if (d == 'a')
	{
		new.model = update.alpha(model);
	}else if (d == 'p')
	{
		new.model = update.p(model,allele.freq);
	}else if (d == 'n')
	{
		new.model = update.nu(model);
	}else
	{
		print("No draw!");	
		exit();
	}

	return(new.model);	
}
	

### metropolis hastings

mh.mcmc = function(num.comp, allele.freq, data.set,num.iter = 1000,thin=10)
{
	current  = init.model(num.comp,allele.freq)
	curr.llk = calc.llk(data.set,current,allele.freq)
	curr.lpr = calc.prior(current);

	out = rep(list(),num.iter/thin)


	cat(paste("Component: ",num.comp,"\n",sep=""));
	pb <- txtProgressBar(1,num.iter)
	inc = round(num.iter/20)
	for (i in 1:num.iter)	
	{
		proposed = current;
		proposed = update.model.mcmc(current)

		prop.llk = calc.llk(data.set,proposed,allele.freq)
		prop.lpr = calc.prior(proposed);
		
		log.llk  = prop.llk - curr.llk;
		log.pri  = prop.lpr - curr.lpr;
		log.prop = calc.proposal(proposed,current)

		alpha = log.llk + log.pri + log.prop;

		u = runif(1);

		if (log(u) < min(alpha,0))
		{
			current = proposed;
			curr.llk = prop.llk;
			curr.lpr = prop.lpr;
		}
		if (i %% thin ==0)
		{
		#	print(c(i,curr.llk+curr.lpr))
			out[[i/thin]] = list(current,curr.llk,curr.lpr);
		}
		setTxtProgressBar(pb, i)
	}
	return(out)
}

dbeta.binom.zi = function(K,N,p,nu)
{
	min.val = 0.01;
	min.int = 2;

	nu.prime = nu+max(1/p,1/(1-p));

	alpha = p*nu.prime;
	beta  = (1-p)*nu.prime

	
	prob = dbetabinom.ab(K, N, alpha, beta);
	return(prob);
}

### need to use this to class(xx) <- "ps"

`[.ps` <- function(x, i) {
  class(x) <- "list"
  structure(x[i], class="ps")
}

`>.ps` <- function(e1, e2) {
  length(e1[[1]]) > length(e2[[1]])
}

`==.ps` <- function(e1, e2) {
  length(e1[[1]]) == length(e2[[1]])
}


