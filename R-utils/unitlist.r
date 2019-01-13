#==========================================================================================#
#==========================================================================================#
#     This list contains the list of all commonly used units written in mathematical       #
# format.                                                                                  #
#------------------------------------------------------------------------------------------#
untab <<- list( cm          = "c*m"
              , cm2m        = "c*m^2*m"
              , cm2om2      = "c*m^2*m^{-2}"
              , cm2opl      = "c*m^2*p*l*a*n*t^{-1}"
              , cm2om2oyr   = "c*m^2*m^{-2}*y*r^{-1}"
              , day         = "d*a*y"
              , deg         = "degree"
              , degC        = "degree*C"
              , degE        = "degree*E"
              , degN        = "degree*N"
              , degS        = "degree*S"
              , degW        = "degree*W"
              , empty       = "phantom(1)-phantom(1)"
              , gcokgw      = "g[C]^phantom(1)*k*g[W]^{-1}"
              , gcokgcbio   = "g[C]^phantom(1)*k*g[C[b*i*o]]^{-1}"
              , Gjom2       = "G*J^phantom(1)*m^{-2}"
              , gmt         = "G*M*T"
              , gom3        = "g^phantom(1)*m^{-3}"
              , gocm3       = "g^phantom(1)*c*m^{-3}"
              , gwokg       = "g[W]^phantom(1)*k*g^{-1}"
              , gwokgodec   = "g[W]^phantom(1)*k*g^{-1}*d*e*c*a*d*e^{-1}"
              , hpa         = "h*P*a"
              , hpaodec     = "h*P*a^phantom(1)*d*e*c*a*d*e^{-1}"
              , hr          = "h*o*u*r"
              , jokgom2     = "J^phantom(1)*k*g^{-1}*m^{-2}"
              , jokgopl     = "J^phantom(1)*k*g^{-1}*p*l^{-1}"
              , jom2        = "J^phantom(1)*m^{-2}"
              , K           = "K"
              , Kodec       = "K^phantom(1)*d*e*c*a*d*e^{-1}"
              , k           = "K"
              , kg          = "k*g"
              , kgcokgc     = "k*g[C]^phantom(1)*k*g[C]^{-1}"
              , kgcokgn     = "k*g[C]^phantom(1)*k*g[N]^{-1}"
              , kgcokgp     = "k*g[C]^phantom(1)*k*g[P]^{-1}"
              , kgcom2      = "k*g[C]^phantom(1)*m^{-2}"
              , kgcom2l     = "k*g[C]^phantom(1)*m[l*e*a*f]^{-2}"
              , kgcom2loyr  = "k*g[C]^phantom(1)*m[l*e*a*f]^{-2}*y*r^{-1}"
              , kgcom2oday  = "k*g[C]^phantom(1)*m^{-2}*d*a*y^{-1}"
              , kgcom2oyr   = "k*g[C]^phantom(1)*m^{-2}*y*r^{-1}"
              , kgcom3      = "k*g[C]^phantom(1)*m^{-3}"
              , kgcopl      = "k*g[C]^phantom(1)*p*l*a*n*t^{-1}"
              , kgcoployr   = "k*g[C]^phantom(1)*p*l*a*n*t^{-1}*y*r^{-1}"
              , kgokg       = "k*g^phantom(1)*k*g^{-1}"
              , kgom2       = "k*g^phantom(1)*m^{-2}"
              , kgom3       = "k*g^phantom(1)*m^{-3}"
              , kgnokg      = "k*g[N]^phantom(1)*k*g^{-1}"
              , kgnom2      = "k*g[N]^phantom(1)*m^{-2}"
              , kgpokg      = "k*g[P]^phantom(1)*k*g^{-1}"
              , kgwokg      = "k*g[W]^phantom(1)*k*g^{-1}"
              , kgpom2      = "k*g[P]^phantom(1)*m^{-2}"
              , kgwom2      = "k*g[W]^phantom(1)*m^{-2}"
              , kgwom2l     = "k*g[W]^phantom(1)*m[l*e*a*f]^{-2}"
              , kgwom2oday  = "k*g[W]^phantom(1)*m^{-2}*d*a*y^{-1}"
              , kgwom2ohr   = "k*g[W]^phantom(1)*m^{-2}*h*r^{-1}"
              , kgwom2os    = "k*g[W]^phantom(1)*m^{-2}*s^{-1}"
              , kgwom2loday = "k*g[W]^phantom(1)*m[l*e*a*f]^{-2}*d*a*y^{-1}"
              , kgwom3oday  = "k*g[W]^phantom(1)*m^{-3}*d*a*y^{-1}"
              , kgwoploday  = "k*g[W]^phantom(1)*p*l*a*n*t^{-1}*d*a*y^{-1}"
              , km          = "k*m"
              , Kmos        = "K^phantom(1)*m^phantom(1)*s^{-1}"
              , m           = "m"
              , m2          = "m^2"
              , mm          = "m*m"
              , mmoyrodec   = "m*m^phantom(1)*y*r^{-1}*d*e*c*a*d*e^{-1}"
              , mm2okgw     = "m*m^2*k*g[W]^{-1}"
              , mmoday      = "m*m^phantom(1)*d*a*y^{-1}"
              , mmolom2os   = "m*m*o*l^phantom(1)*m^{-2}*s^{-1}"
              , mmomo       = "m*m^phantom(1)*m*o*n*t*h^{-1}"
              , molom2      = "m*o*l^phantom(1)*m^{-2}"
              , molom2l     = "m*o*l^phantom(1)*m[l*e*a*f]^{-2}"
              , molcom2     = "m*o*l[C]^phantom(1)*m^{-2}"
              , month       = "m*o*n*t*h"
              , mos         = "m^phantom(1)*s^{-1}"
              , mosodec     = "m^phantom(1)*s^{-1}*d*e*c*a*d*e^{-1}"
              , m2om2       = "m^2*m^{-2}"
              , m2om3       = "m^2*m^{-3}"
              , m2opl       = "m^2*p*l^{-1}"
              , m2lokgc     = "m[l*e*a*f]^2*k*g[C]^{-1}"
              , m2lom2      = "m[l*e*a*f]^2*m^{-2}"
              , m2lopl      = "m[l*e*a*f]^2*p*l^{-1}"
              , m2pom2      = "m[l*e*a*f+w*o*o*d]^2*m^{-2}"
              , m2popl      = "m[l*e*a*f+w*o*o*d]^2*p*l^{-1}"
              , m2wom2      = "m[w*o*o*d]^2*m^{-2}"
              , m2wopl      = "m[w*o*o*d]^2*p*l^{-1}"
              , m3oha       = "m^3*h*a^{-1}"
              , m3om2       = "m^3*m^{-2}"
              , m3wom3      = "m[W]^3*m^{-3}"
              , Mgwom2      = "M*g[W]^phantom(1)*m^{-2}"
              , Mjom2       = "M*J^phantom(1)*m^{-2}"
              , mmoyr       = "m*m^phantom(1)*y*r^{-1}"
              , mpa         = "M*P*a"
              , nmo.090     = "m*o*n*t*h*s*phantom(1)*\"|\"*phantom(1)*bar(R) < 90*m*m^phantom(1)*m*o^{-1}"
              , nmo.100     = "m*o*n*t*h*s*phantom(1)*\"|\"*phantom(1)*bar(R) < 100*m*m^phantom(1)*m*o^{-1}"
              , nmo.110     = "m*o*n*t*h*s*phantom(1)*\"|\"*phantom(1)*bar(R) < 110*m*m^phantom(1)*m*o^{-1}"
              , nmo.120     = "m*o*n*t*h*s*phantom(1)*\"|\"*phantom(1)*bar(R) < 120*m*m^phantom(1)*m*o^{-1}"
              , nmo.wdef    = "m*o*n*t*h*s*phantom(1)*\"|\"*phantom(1)*bar(W*D) > 10*m*m^phantom(1)*m*o^{-1}"
              , nmolomol    = "n*m*o*l^phantom(1)*m*o*l^{-1}"
              , oneom       = "m^{-1}"
              , oneoha      = "h*a^{-1}"
              , oneom2      = "m^{-2}"
              , oneoyr      = "y*r^{-1}"
              , pa          = "P*a"
              , pc          = "'%'"
              , pcbio       = "'%'[b*i*o]"
              , pcagboyr    = "'%'[A*G*B]^phantom(1)*y*r^{-1}"
              , pcagbo50yr  = "'%'[A*G*B]^phantom(1)*group(\"(\",50*y*r,\")\")^{-1}"
              , pcbaoyr     = "'%'[B*A]^phantom(1)*y*r^{-1}"
              , pcbiooyr    = "'%'[b*i*o]^phantom(1)*y*r^{-1}"
              , pcdbhoyr    = "'%'[D*B*H]^phantom(1)*y*r^{-1}"
              , pcetoyr     = "'%'[E*T]^phantom(1)*y*r^{-1}"
              , pceto50yr   = "'%'[E*T]^phantom(1)*group(\"(\",50*y*r,\")\")^{-1}"
              , pcoyr       = "'%'^phantom(1)*y*r^{-1}"
              , pcpopoyr    = "'%'[p*o*p]^phantom(1)*y*r^{-1}"
              , pcsat       = "'%'[S*a*t]"
              , plom2       = "p*l*a*n*t^phantom(1)*m^{-2}"
              , s           = "s"
              , thkm2oyr    = "10^{3*phantom(1)}*k*m^{2}*y*r^{-1}"
              , umolcom2os  = "mu*m*o*l[C]^phantom(1)*m^{-2}*s^{-1}"
              , umolokgcos  = "mu*m*o*l^phantom(1)*k*g[C]^{-1}*s^{-1}"
              , umolom2     = "mu*m*o*l^phantom(1)*m^{-2}"
              , umolom2os   = "mu*m*o*l^phantom(1)*m^{-2}*s^{-1}"
              , umolom2l    = "mu*m*o*l^phantom(1)*m[l*e*a*f]^{-2}"
              , umolom2los  = "mu*m*o*l^phantom(1)*m[l*e*a*f]^{-2}*s^{-1}"
              , umolcom2    = "mu*m*o*l[C]^phantom(1)*m^{-2}"
              , umolcomol   = "mu*m*o*l[C]^phantom(1)*m*o*l^{-1}"
              , umolomol    = "mu*m*o*l^phantom(1)*m*o*l^{-1}"
              , utc         = "U*T*C"
              , wom2        = "W^phantom(1)*m^{-2}"
              , wom2l       = "W^phantom(1)*m[l*e*a*f]^{-2}"
              , wom2odec    = "W^phantom(1)*m^{-2}*d*e*c*a*d*e^{-1}"
              , wopl        = "W^phantom(1)*p*l*a*n*t^{-1}"
              , yr          = "y*r"
              )#end list
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This list contains the inverse of all commonly used units written in mathematical    #
# format.  The inverse is useful for probability distribution functions (units of the      #
# density function are the inverse of the actual units of the quantity, so the integral    #
# becomes dimensionless.                                                                   #
#------------------------------------------------------------------------------------------#
i.untab <<- list( cm          = "c*m^{-1}"
                , cm2m        = "c*m^{-2}*m^{-1}"
                , cm2om2      = "m^2*c*m^{-2}"
                , cm2opl      = "p*l*a*n*t^phantom(1)*c*m^{-2}"
                , cm2om2oyr   = "m^2*y*r^phantom(1)*c*m^{-2}"
                , day         = "d*a*y^{-1}"
                , deg         = "degree^{-1}"
                , degC        = "degree*C^{-1}"
                , degE        = "degree*E^{-1}"
                , degN        = "degree*N^{-1}"
                , degS        = "degree*S^{-1}"
                , degW        = "degree*W^{-1}"
                , empty       = "phantom(1)-phantom(1)"
                , gcokgw      = "k*g[W]^phantom(1)*g[C]^{-1}"
                , gcokgcbio   = "k*g[C[b*i*o]]^phantom(1)*g[C]^{-1}"
                , gmt         = "G*M*T^{-1}"
                , gom3        = "m^3*g^{-1}"
                , gocm3       = "c*m^3*g^{-1}"
                , gwokg       = "k*g^phantom(1)*g[W]^{-1}"
                , gwokgodec   = "k*g^phantom(1)*d*e*c*a*d*e^phantom(1)*k*g[W]^{-1}"
                , hpa         = "h*P*a^{-1}"
                , hpaodec     = "d*e*c*a*d*e^phantom(1)*h*P*a^{-1}"
                , hr          = "h*o*u*r^{-1}"
                , jokgom2     = "k*g^phantom(1)*m^{-2}*J^{-1}"
                , jokgopl     = "k*g^phantom(1)*pl^phantom(1)*J^{-1}"
                , jom2        = "m^{2}*J^{-1}"
                , K           = "K^{-1}"
                , Kodec       = "d*e*c*a*d*e^phantom(1)*K^{-1}"
                , k           = "K^{-1}"
                , kg          = "k*g^{-1}"
                , kgcokgc     = "k*g[C]^phantom(1)*k*g[C]^{-1}"
                , kgcokgn     = "k*g[N]^phantom(1)*k*g[C]^{-1}"
                , kgcokgp     = "k*g[P]^phantom(1)*k*g[C]^{-1}"
                , kgcom2      = "m^{2}*k*g[C]^{-1}"
                , kgcom2l     = "m[l*e*a*f]^{2}*k*g[C]^{-1}"
                , kgcom2loyr  = "m[l*e*a*f]^{2}*y*r^phantom(1)*k*g[C]^{-1}"
                , kgwom2oday  = "m^{2}*d*a*y^phantom(1)*k*g[C]^{-1}"
                , kgcom2oyr   = "m^{2}*y*r^phantom(1)*k*g[C]^{-1}"
                , kgcom3      = "m^{3}*k*g[C]^{-1}"
                , kgcopl      = "p*l*a*n*t^phantom(1)*k*g[C]^{-1}"
                , kgcoployr   = "p*l*a*n*t^phantom(1)*y*r^phantom(1)*k*g[C]^{-1}"
                , kgokg       = "k*g^phantom(1)*k*g^{-1}"
                , kgom2       = "m^{2}*k*g^{-1}"
                , kgom3       = "m^{3}*k*g^{-1}"
                , kgnokg      = "k*g^phantom(1)*k*g[N]^{-1}"
                , kgnom2      = "m^{2}*k*g[N]^{-1}"
                , kgpokg      = "k*g^phantom(1)*k*g[P]^{-1}"
                , kgwokg      = "k*g^phantom(1)*k*g[W]^{-1}"
                , kgpom2      = "m^{2}*k*g[P]^{-1}"
                , kgwom2      = "m^{2}*k*g[W]^{-1}"
                , kgwom2l     = "m[l*e*a*f]^{2}*k*g[W]^{-1}"
                , kgwom2oday  = "m^{2}*d*a*y^phantom(1)*k*g[W]^{-1}"
                , kgwom2ohr   = "m^{2}*h*r^phantom(1)*k*g[W]^{-1}"
                , kgwom2os    = "m^{2}*s^phantom(1)*k*g[W]^{-1}"
                , kgwom2loday = "m[l*e*a*f]^{2}*d*a*y^phantom(1)*k*g[W]^{-1}"
                , kgwom3oday  = "m^{-3}*d*a*y^phantom(1)*k*g[W]^{-1}"
                , kgwoploday  = "p*l*a*n*t^phantom(1)*d*a*y^phantom(1)*k*g[W]^{-1}"
                , km          = "k*m^{-1}"
                , Kmos        = "s^phantom(1)*K^{-1}*m^{-1}"
                , m           = "m^{-1}"
                , m2          = "m^{-2}"
                , mm          = "m*m^{-1}"
                , mmoyrodec   = "y*r^phantom(1)*d*e*c*a*d*e^phantom(1)*m*m^^{-1}"
                , mm2okgw     = "k*g[W]^phantom(1)*m*m^{-2}"
                , mmoday      = "d*a*y^phantom(1)*m*m^{-1}"
                , mmolom2os   = "m^{2}*s^phantom(1)*m*m*o*l^{-1}"
                , mmomo       = "m*o*n*t*h^phantom(1)*^m*m{-1}"
                , molom2      = "m^{2}*m*o*l^{-1}"
                , molom2l     = "m[l*e*a*f]^{2}*m*o*l^{-1}"
                , molcom2     = "m^{2}*m*o*l[C]^{-1}"
                , month       = "m*o*n*t*h^{-1}"
                , mos         = "s^phantom(1)*m^{-1}"
                , mosodec     = "s^phantom(1)*d*e*c*a*d*e^phantom(1)*m^{-1}"
                , m2om2       = "m^2*m^{-2}"
                , m2om3       = "m^3*m^{-2}"
                , m2opl       = "p*l^phantom(1)*m^{-2}"
                , m2lokgc     = "k*g[C]^phantom(1)*m[l*e*a*f]^{-2}"
                , m2lom2      = "m^2*m[l*e*a*f]^{-2}"
                , m2lopl      = "p*l^phantom(1)*m[l*e*a*f]^{-2}"
                , m2pom2      = "m^2*m[l*e*a*f+w*o*o*d]^{-2}"
                , m2popl      = "p*l^phantom(1)*m[l*e*a*f+w*o*o*d]^{-2}"
                , m2wom2      = "m^2*m[w*o*o*d]^{-2}"
                , m2wopl      = "p*l^phantom(1)*m[w*o*o*d]^{-2}"
                , m3oha       = "h*a^phantom(1)*m^{-3}"
                , m3om2       = "m^2*m^{-3}"
                , m3wom3      = "m^3*m[W]^{-3}"
                , Mjom2       = "m^{2}*M*J^{-1}"
                , mmoyr       = "y*r^phantom(1)*m*m^{-1}"
                , mpa         = "M*P*a^{-1}"
                , nmo.090     = "m*o*n*t*h*s*{-1}*\"|\"*phantom(1)*bar(R) < 90*m*m^phantom(1)*m*o^{-1}"
                , nmo.100     = "m*o*n*t*h*s*{-1}*\"|\"*phantom(1)*bar(R) < 100*m*m^phantom(1)*m*o^{-1}"
                , nmo.110     = "m*o*n*t*h*s*{-1}*\"|\"*phantom(1)*bar(R) < 110*m*m^phantom(1)*m*o^{-1}"
                , nmo.120     = "m*o*n*t*h*s*{-1}*\"|\"*phantom(1)*bar(R) < 120*m*m^phantom(1)*m*o^{-1}"
                , nmo.wdef    = "m*o*n*t*h*s*{-1}*\"|\"*phantom(1)*bar(W*D) > 10*m*m^phantom(1)*m*o^{-1}"
                , nmolomol    = "m*o*l^phantom(1)*n*m*o*l^{-1}"
                , oneom       = "m"
                , oneoha      = "h*a"
                , oneom2      = "m^{2}"
                , oneoyr      = "y*r"
                , pa          = "P*a^{-1}"
                , pc          = "'%'^{-1}"
                , pcbio       = "'%'[b*i*o]^{-1}"
                , pcagboyr    = "y*r^phantom(1)*'%'[A*G*B]^^{-1}"
                , pcagbo50yr  = "group(\"(\",50*y*r,\")\")^phantom(1)*'%'[A*G*B]^{-1}"
                , pcbaoyr     = "y*r^phantom(1)*'%'[B*A]^{-1}"
                , pcbiooyr    = "y*r^phantom(1)*'%'[b*i*o]^{-1}"
                , pcdbhoyr    = "y*r^phantom(1)*'%'[D*B*H]^{-1}"
                , pcetoyr     = "y*r^phantom(1)*'%'[E*T]^{-1}"
                , pceto50yr   = "group(\"(\",50*y*r,\")\")^phantom(1)*'%'[E*T]^{-1}"
                , pcoyr       = "y*r^phantom(1)*'%'^{-1}"
                , pcpopoyr    = "y*r^phantom(1)*'%'[p*o*p]^{-1}"
                , pcsat       = "'%'[S*a*t]^{-1}"
                , plom2       = "m^2*p*l*a*n*t^{-1}"
                , s           = "s^{-1}"
                , thkm2oyr    = "y*r^phantom(1)*10^{-3*phantom(1)}*k*m^{-2}"
                , umolcom2os  = "m^2*s^phantom(1)*mu*m*o*l[C]^{-1}"
                , umolokgcos  = "k*g[C]^phantom(1)*s^phantom(1)*mu*m*o*l^{-1}"
                , umolom2     = "m^2*mu*m*o*l^{-1}"
                , umolom2os   = "m^2*s^phantom(1)*mu*m*o*l^{-1}"
                , umolom2l    = "m[l*e*a*f]^2*mu*m*o*l^{-1}"
                , umolom2los  = "m[l*e*a*f]^2*s^phantom(1)*mu*m*o*l^{-1}"
                , umolcom2    = "m^2*mu*m*o*l[C]^{-1}"
                , umolcomol   = "m*o*l^phantom(1)*mu*m*o*l[C]^{-1}"
                , umolomol    = "m*o*l^phantom(1)*mu*m*o*l^{-1}"
                , utc         = "U*T*C^{-1}"
                , wom2        = "m^2*W^{-1}"
                , wom2l       = "m[l*e*a*f]^{2}*W^phantom(1)"
                , wom2odec    = "m^{2}*d*e*c*a*d*e^phantom(1)*W^{-1}"
                , wopl        = "p*l*a*n*t^phantom(1)*W^{-1}"
                , yr          = "y*r^{-1}"
                )#end list
#==========================================================================================#
#==========================================================================================#
