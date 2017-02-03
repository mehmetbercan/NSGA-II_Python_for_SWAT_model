# SWAT MAKEFILE FOR INTEL COMPILER on LINUX
#
F95=ifort

FFLAGS= -O1 -ipo -debug none -fpe:0 -fp-model precise -zero
#FFLAGS= -fast -ipo -fpe:0

LDFLAGS= -L/usr/local/lib $(FFLAGS) -static -static-intel
#LDFLAGS= -L/usr/local/lib $(FFLAGS)

MOD=modparm.o

EXE=swat2012_627


SRC = addh.f albedo.f allocate_parms.f alph.f anfert.f apex_day.f apply.f ascrv.f atri.f aunif.f autoirr.f aveval.f bacteria.f biofilm.f biozone.f bmp_det_pond.f bmpfixed.f bmpinit.f bmp_ri_pond.f bmp_sand_filter.f bmp_sed_pond.f bmp_wet_pond.f buffer.f burnop.f canopyint.f caps.f carbon_new.f carbon_zhang2.f90 cfactor.f chkcst.f clgen.f clicon.f command.f conapply.f confert.f crackflow.f crackvol.f curno.f dailycn.f decay.f depstor.f distrib_bmps.f dormant.f drains.f dstn1.f ee.f eiusle.f enrsb.f erfc.f estimate_ksat.f etact.f etpot.f expo.f fert.f filter.f filtw.f finalbal.f gcycl.f getallo.f grass_wway.f graze.f grow.f gwmod_deep.f gwmod.f gw_no3.f gwnutr.f h2omgt_init.f harvestop.f harvgrainop.f harvkillop.f header.f headout.f hhnoqual.f hhwatqual.f hmeas.f hruaa.f hruallo.f hruday.f90 hrumon.f hrupond.f hrupondhr.f hruyr.f hydroinit.f icl.f impndaa.f impndday.f impnd_init.f impndmon.f impndyr.f irrigate.f irr_rch.f irr_res.f irrsub.f jdt.f killop.f lakeq.f latsed.f layersplit.f lwqdef.f main.f NCsed_leach.f90 ndenit.f newtillmix.f nfix.f nitvol.f nlch.f nminrl.f noqual.f npup.f nrain.f nup.f nuts.f openwth.f operatn.f orgncswat.f orgn.f origtile.f ovr_sed.f percmacro.f percmain.f percmicro.f pestlch.f pestw.f pesty.f pgen.f pgenhr.f pkq.f plantmod.f plantop.f pmeas.f pminrl2.f pminrl.f pond.f pondhr.f pothole.f potholehr.f print_hyd.f psed.f qman.f ran1.f rchaa.f rchday.f rchinit.f rchmon.f rchuse.f rchyr.f reachout.f readatmodep.f readbsn.f readchm.f readcnst.f readfcst.f readfert.f readfig.f readfile.f readgw.f readhru.f readinpt.f readlup.f readlwq.f readmgt.f readmon.f readops.f readpest.f readplant.f readpnd.f readres.f readrte.f readru.f readsdr.f readsepticbz.f readseptwq.f readsno.f readsol.f readsub.f readswq.f readtill.f readurban.f readwgn.f readwus.f readwwq.f readyr.f reccnst.f recday.f rechour.f recmon.f recyear.f regres.f resbact.f resetlu.f res.f reshr.f resinit.f resnut.f rewind_init.f rhgen.f rootfr.f route.f routels.f routeunit.f routres.f rsedaa.f rseday.f rsedmon.f rsedyr.f rtbact.f rtday.f rteinit.f rthmusk.f rthpest.f rthr.f rthsed.f rtmusk.f rtout.f rtpest.f rtsed_bagnold.f rtsed.f rtsed_kodatie.f rtsed_Molinas_Wu.f rtsed_yangsand.f sat_excess.f saveconc.f save.f sched_mgt.f schedule_ops.f sim_initday.f sim_inityr.f simulate.f slrgen.f smeas.f snom.f soil_chem.f soil_par.f soil_phys.f soil_write.f solp.f solt.f std1.f std2.f std3.f stdaa.f storeinitial.f structure.f subaa.f subbasin.f subday.f submon.f substor.f sub_subbasin.f subwq.f subyr.f sumhyd.f sumv.f surface.f surfst_h2o.f surfstor.f surq_daycn.f surq_greenampt.f swbl.f sweep.f swu.f tair.f tgen.f theta.f tillfactor.f tillmix.f tmeas.f tran.f transfer.f tstr.f ttcoef.f ttcoef_wway.f urban.f urbanhr.f urb_bmp.f varinit.f vbl.f virtual.f volq.f washp.f watbal.f water_hru.f watqual2.f watqual.f wattable.f watuse.f weatgn.f wetlan.f wmeas.f wndgen.f writeaa.f writea.f writed.f writem.f xmon.f ysed.f zero0.f zero1.f zero2.f zeroini.f zero_urbn.f

OBJ=$(SRC:.f=.o)

all: $(EXE)

clean:
	rm -f *.o $(MOD) $(EXE)

$(MOD): modparm.f
	$(F95) -c $(FFLAGS) $<

$(OBJ): $(MOD)

$(EXE): $(OBJ) modparm.o
	$(F95) $(LDFLAGS) $(OBJ) -o $@

.f.o: modparm.f
	$(F95) -c $(FFLAGS) $< -o $@

