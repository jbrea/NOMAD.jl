module NOMAD
using Cxx, MathProgBase

const depfile = joinpath(dirname(@__FILE__),"..","deps","deps.jl")
if isfile(depfile)
    include(depfile)
else
    error("NOMAD not properly installed. Please run Pkg.build(\"NOMAD\")")
end
const nomaddir = joinpath(dirname(@__FILE__),"..","deps", "src", "nomad")
addHeaderDir(joinpath(nomaddir, "src"), kind = C_System)
addHeaderDir(joinpath(nomaddir, "ext", "sgtelib", "src"), kind = C_System)
Libdl.dlopen(libnomad, Libdl.RTLD_GLOBAL)
Libdl.dlopen(libsgtelib, Libdl.RTLD_GLOBAL)

cxx"""
	#include "nomad.hpp"
	using namespace std;
"""

import MathProgBase.SolverInterface
import MathProgBase.SolverInterface.optimize!

struct Evaluator
	ndims::Int
	cppclassname::String
	callbackfunction::Function
	callbackconstraints::Array{Function, 1}
end
mutable struct CxxID
	ev::Int
	main::Int
end
const CXXID = CxxID(0, 0)
function Evaluator(func, constraints, ndims)
	CXXID.ev += 1
	cbfuncstring = "function evalfunc_$(CXXID.ev)($(join(["_x$i" for i in 1:ndims], ", ")))\n"
	for i in 1:ndims
		cbfuncstring *= "$(Symbol(:x, i)) = unsafe_load($(Symbol(:_x, i)))\n"
	end
 	cbfuncstring *= "Float64(Main.$(func)($(join(["x$i" for i in 1:ndims], ", "))))\nend"
# 	println(cbfuncstring)
	cbfunc = eval(parse(cbfuncstring))
	cppclassname = "NOMAD_Evaluator_$(CXXID.ev)"
	cxxstring = """
		class $cppclassname : public NOMAD::Evaluator {
			public:
				$cppclassname  ( const NOMAD::Parameters & p ) : NOMAD::Evaluator ( p ) {}
		~$cppclassname ( void ) {}
		bool eval_x ( NOMAD::Eval_Point   & x          ,
					  const NOMAD::Double & h_max      ,
					  bool                & count_eval  ) const
			{
			x.set_bb_output  ( 0 , \$($(cbfunc))($(join(["x[$i].value()" 
													 for i in 0:ndims-1],
													","))));
			count_eval = true; 
			return true;
			}        
		};"""
#  	println(cxxstring)
	eval(Cxx.process_cxx_string(cxxstring))
	Evaluator(ndims, cppclassname, cbfunc, Function[])
end

struct Opt
	evaluator::Evaluator
	cppfunc::Function
end
function Opt(ev; kargs...)
	paramstring = parameterstring(ev.ndims, kargs)
	CXXID.main += 1
	funcname = Symbol(:main_, CXXID.main)
	cxxstring = """int $funcname() {
		NOMAD::Display out ( std::cout );
        out.precision ( NOMAD::DISPLAY_PRECISION_STD );
        try {
			NOMAD::Parameters p (out); 
			p.set_DIMENSION ($(ev.ndims));
			vector<NOMAD::bb_output_type> bbot (1);
			bbot[0] = NOMAD::OBJ;
			p.set_BB_OUTPUT_TYPE( bbot ); $paramstring; p.check();
			$(ev.cppclassname) ev (p);
			NOMAD::Mads mads (p, &ev);
			mads.run();
		}
		catch ( exception & e ) {
           cerr << "\\nNOMAD has been interrupted (" << e.what() << ")\\n\\n\";
        }
        NOMAD::Slave::stop_slaves ( out );
        NOMAD::end();
        return EXIT_SUCCESS;
       }
	"""
# 	println(cxxstring)
	eval(Cxx.process_cxx_string(cxxstring))
	cppfunc() = @cxx $funcname()
	Opt(ev, cppfunc)
end

function preprocessarg(ndims, key, arg)
	if key == :UPPER_BOUND || key == :LOWER_BOUND || key == :X0
		s = "NOMAD::Point point_$key ( $ndims );"
		for i in 0:ndims - 1
			s *= "point_$key[$i] = $(arg[i+1]);"
		end
		s *= "p.set_$key( point_$key );"
		s
	else
		"p.set_$key( $arg );"
	end
end
function parameterstring(ndims, kargs)
	s = ""
	for (key, arg) in kargs
		s *= preprocessarg(ndims, key, arg)
	end
	s
end

function help(key = "")
	bindir = joinpath(nomaddir, "bin", "nomad")
	s = readstring(`$bindir --help $key`)
	print(replace(s, "\$NOMAD_HOME", normpath(nomaddir)))
end


struct NOMADSolver <: SolverInterface.AbstractMathProgSolver
	
end

mutable struct NOMADMathProgModel <: SolverInterface.AbstractNonlinearModel

end

end # module
