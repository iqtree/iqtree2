add_library(model
modelmarkov.cpp modelmarkov.h
modelbin.cpp modelbin.h
modeldna.cpp modeldna.h
modeldnaerror.cpp modeldnaerror.h
modelfactory.cpp modelfactory.h
modelprotein.cpp modelprotein.h
modelset.cpp modelset.h
modelsubst.cpp modelsubst.h
partitionmodel.cpp partitionmodel.h
partitionmodelplen.cpp partitionmodelplen.h
rategamma.cpp rategamma.h
rategammainvar.cpp rategammainvar.h
ratecontinuousgamma.cpp ratecontinuousgamma.h
ratecontinuousgammainvar.cpp ratecontinuousgammainvar.h
rateheterogeneity.cpp rateheterogeneity.h
rateinvar.cpp rateinvar.h
ratemeyerdiscrete.cpp ratemeyerdiscrete.h
ratemeyerhaeseler.cpp ratemeyerhaeseler.h
ratekategory.cpp ratekategory.h
ratefree.cpp ratefree.h
ratefreeinvar.cpp ratefreeinvar.h
rateheterotachy.cpp rateheterotachy.h
rateheterotachyinvar.cpp rateheterotachyinvar.h
modelcodon.cpp modelcodon.h
modelmorphology.cpp modelmorphology.h
modelmixture.cpp modelmixture.h
modelunrest.cpp modelunrest.h
modelliemarkov.cpp modelliemarkov.h
modelpomo.cpp modelpomo.h
modelpomomixture.cpp modelpomomixture.h
modelfactorymixlen.cpp modelfactorymixlen.h
)

target_link_libraries(model utils)
