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

function createcbfunc(id, func, ndims)
    cbfuncstring = "function evalfunc_$id($(join(["_x$i" for i in 1:ndims], ", ")))\n"
    for i in 1:ndims
        cbfuncstring *= "$(Symbol(:x, i)) = unsafe_load($(Symbol(:_x, i)))\n"
    end
    cbfuncstring *= "Float64(Main.$(func)($(join(["x$i" for i in 1:ndims], ", "))))\nend"
    cbfunc = eval(parse(cbfuncstring))
end

struct Evaluator
    ndims::Int
    cppclassname::String
    callbackfunction::Function
    callbackconstraints::Array{Tuple{Symbol, Function}, 1}
end
mutable struct CxxID
    ev::Int
    main::Int
end
const CXXID = CxxID(0, 0)
function Evaluator(func, constraints, ndims; printcxx = false)
    CXXID.ev += 1
    cppclassname = "NOMAD_Evaluator_$(CXXID.ev)"
    objfunc = createcbfunc(CXXID.ev, func, ndims)
    argstring = join(["x[$i].value()" for i in 0:ndims-1], ",")
    constrfuncs = Tuple{Symbol, Function}[]
    constrstring = ""
    for (i, constr) in enumerate(constraints)
        if typeof(constr) <: Function
            constr = (:EB, constr)
        end
        push!(constrfuncs, (constr[1], 
                            createcbfunc("$(CXXID.ev)_c$i", constr[2], ndims)))
        constrstring *= "\n    x.set_bb_output  ( $i , \$($(constrfuncs[i][2]))($argstring));"
    end
    cxxstring = """
        class $cppclassname : public NOMAD::Evaluator {
            public:
                $cppclassname  ( const NOMAD::Parameters & p ) : NOMAD::Evaluator ( p ) {}
        ~$cppclassname ( void ) {}
        bool eval_x ( NOMAD::Eval_Point   & x          ,
                      const NOMAD::Double & h_max      ,
                      bool                & count_eval  ) const
            {
            x.set_bb_output  ( 0 , \$($(objfunc))($argstring)); $constrstring;
            count_eval = true; 
            return true;
            }        
        };"""
    if printcxx; println(cxxstring); end
    eval(Cxx.process_cxx_string(cxxstring))
    Evaluator(ndims, cppclassname, objfunc, constrfuncs)
end
export Evaluator

struct Opt
    evaluator::Evaluator
    cppfunc::Function
end
function Opt(ev; printcxx = false, kargs...)
    paramstring = parameterstring(ev.ndims, kargs)
    CXXID.main += 1
    funcname = Symbol(:main_, CXXID.main)
    bboutdim = length(ev.callbackconstraints) + 1
    bboutdefstring = "bbot[0] = NOMAD::OBJ;"
    for (i, (typ, f)) in enumerate(ev.callbackconstraints)
        bboutdefstring *= "\n        bbot[$i] = NOMAD::$typ;"
    end
    cxxstring = """int $funcname(double* bf_x, double* bf_f, 
                                 NOMAD::Stats *stats) {
        NOMAD::Display out ( std::cout );
        out.precision ( NOMAD::DISPLAY_PRECISION_STD );
        try {
            NOMAD::Parameters p (out); 
            p.set_DIMENSION ($(ev.ndims));
            vector<NOMAD::bb_output_type> bbot ($bboutdim);
            $bboutdefstring;
            p.set_BB_OUTPUT_TYPE( bbot ); $paramstring; 
            p.check();
            $(ev.cppclassname) ev (p);
            NOMAD::Mads mads (p, &ev);
            mads.run();

            const NOMAD::Eval_Point* bf = mads.get_best_feasible();
            bf_f[0] = bf->get_f().value();
            for ( int i = 0; i < $(ev.ndims); i++) {
                bf_x[i] = bf->value(i);
            }
            *stats = mads.get_stats();
        }
        catch ( exception & e ) {
           cerr << "\\nNOMAD has been interrupted (" << e.what() << ")\\n\\n\";
        }
        NOMAD::Slave::stop_slaves ( out );
        NOMAD::end();
        return EXIT_SUCCESS;
       }
    """
    if printcxx; println(cxxstring); end
    eval(Cxx.process_cxx_string(cxxstring))
    cppfunc(bf_x, bf_f, stats) = @cxx $funcname(pointer(bf_x), 
                                                pointer(bf_f),
                                                stats)
    Opt(ev, cppfunc)
end
export Opt

function optimize!(opt)
    bf_x = zeros(opt.evaluator.ndims)
    bf_f = [0.]
    stats = icxx"new NOMAD::Stats();"
    opt.cppfunc(bf_x, bf_f, stats)
    bf_x, bf_f[1], icxx"$stats->get_bb_eval();"
end
export optimize!

function preprocessarg(ndims, key, arg)
    if key == :UPPER_BOUND || key == :LOWER_BOUND || key == :X0
        s = "\n        NOMAD::Point point_$key ( $ndims );"
        for i in 0:ndims - 1
            s *= "\n        point_$key[$i] = $(arg[i+1]);"
        end
        s *= "\n        p.set_$key( point_$key );"
        s
    else
        "\n        p.set_$key( $arg );"
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
    s = "This help message is obtained by calling `nomad --help` (in batch mode).\nNot all features are available in julia.\nPlease open an issue on github, if something doesn't work as expected.\n\n"
    s *= readstring(`$bindir --help $key`)
    print(replace(s, "\$NOMAD_HOME", normpath(nomaddir)))
end


struct NOMADSolver <: SolverInterface.AbstractMathProgSolver
    
end

mutable struct NOMADMathProgModel <: SolverInterface.AbstractNonlinearModel

end

end # module
