###this script checks to see whether error causing VCV matrix to calculate covariance on deeper root affected tanager data
##it does not

load('~/Dropbox/Manuscripts/PLOS Biology Tanager Project/data/stochastic maps for MCC tree/diet.simmaps.RData')

master<-read.csv("~/Dropbox/Manuscripts/PLOS Biology Tanager Project/data/MASTER_tandata.csv")

dietcats<-c("fruit","inv","omn","seed")

for(i in 1:100){

	smap<-diet.simmaps[[i]]
	
	for(j in 1:4){
	
		#trim simmap to preserve branches of interest
		smap.trimmed<-trimSimmap(smap,trim.class=dietcats[j])
		
		#trim tree to tips in subcategory
		tree.trimmed<-drop.tip.simmap(smap,smap$tip.label[which(!smap$tip.label%in%subset(master,diet==dietcats[j])$treename)])
		
		geo.simmap.root=max(nodeHeights(smap.trimmed))
		geo.simmap.trimmed.root=max(nodeHeights(tree.trimmed))

		##roots aren't the same age (as desired); now designing a workaround

		if(round(geo.simmap.root,5)!=round(geo.simmap.trimmed.root,5)){
		eval(parse(text=paste0("print('error with i = ",i,", j = ",j,"')")))
		
		}
		
	}
	
	if(i%%10==0){print(i)}
	
}	

###need to also test habitat and year-round territoriality
load('~/Dropbox/Manuscripts/PLOS Biology Tanager Project/data/stochastic maps for MCC tree/habitat.simmaps.RData')

for(i in 1:100){

	smap<-habitat.simmaps[[i]]
	habcats=c("1","2","3")
	for(j in 1:3){
	
		#trim simmap to preserve branches of interest
		smap.trimmed<-trimSimmap(smap,trim.class=habcats[j])
		
		#trim tree to tips in subcategory
		tree.trimmed<-drop.tip.simmap(smap,smap$tip.label[which(!smap$tip.label%in%subset(master,diet==dietcats[j])$treename)])
		
		geo.simmap.root=max(nodeHeights(smap.trimmed))
		geo.simmap.trimmed.root=max(nodeHeights(tree.trimmed))

		##roots aren't the same age (as desired); now designing a workaround

		if(round(geo.simmap.root,5)!=round(geo.simmap.trimmed.root,5)){
		eval(parse(text=paste0("print('error with i = ",i,", j = ",j,"')")))
		
		}
		
	}
	
	if(i%%10==0){print(i)}
	
}	


load('~/Dropbox/Manuscripts/PLOS Biology Tanager Project/data/stochastic maps for MCC tree/terr.simmaps.RData')

for(i in 1:100){

	smap<-terr.simmaps[[i]]
	
		#trim simmap to preserve branches of interest
		smap.trimmed<-trimSimmap(smap,trim.class="year.round")
		
		#trim tree to tips in subcategory
		tree.trimmed<-drop.tip.simmap(smap,smap$tip.label[which(!smap$tip.label%in%subset(master,diet==dietcats[j])$treename)])
		
		geo.simmap.root=max(nodeHeights(smap.trimmed))
		geo.simmap.trimmed.root=max(nodeHeights(tree.trimmed))

		##roots aren't the same age (as desired); now designing a workaround

		if(round(geo.simmap.root,5)!=round(geo.simmap.trimmed.root,5)){
		eval(parse(text=paste0("print('error with i = ",i,", j = ",j,"')")))
		
		
	}
	
	if(i%%10==0){print(i)}
	
}	
