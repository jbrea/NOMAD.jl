#TODO: 
# 1. return correct status
# 2. is it possible to resume? 


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

function getfuncsigs(funcs, ndims, vectorin = true, resbyref = false)
    args = [Symbol(:x, i) for i in 1:ndims]
    argex = vectorin ? :([$(args...)]) : :($(args...))
    if resbyref
        argex = typeof(argex) <: Tuple ? (argex..., :res) : (argex, :res)
    end
    if typeof(argex) <: Tuple
        [:($f($(argex...))) for f in funcs]
    else
        [:($f($argex)) for f in funcs]
    end
end

function createcbfunc(id, funcs, ndims, vectorin = true, resbyref = false)
    funcsigs = getfuncsigs(funcs, ndims, vectorin, resbyref)
    if resbyref
        funccalls = funcsigs
    else
        mdims = length(funcs)
        funccalls = [:(unsafe_store!(res, $(funcsigs[i]), $i)) for i in 1:mdims]
    end
    cbfuncex = quote
        function $(Symbol(:evalfunc_, id))($([Symbol(:_x, i) 
                                              for i in 1:ndims]...), _res)
            $([:($(Symbol(:x, i)) = unsafe_load($(Symbol(:_x, i))))
               for i in 1:ndims]...)
            res = unsafe_load(_res)
            $(funccalls...)
            nothing
        end
    end
    cbfunc = eval(cbfuncex)
end

struct Evaluator
    ndims::Int
    mdims::Int
    cppclassname::String
    callbackfunction::Function
    constrainttypes::Array{Symbol, 1}
end
mutable struct CxxID
    ev::Int
    main::Int
end
const CXXID = CxxID(0, 0)
function Evaluator(func, constraints, ndims; printcxx = false, vectorin = true)
    mdims = 1 + length(constraints)
    cbfuncs = Function[func]
    constrtypes = Symbol[]
    for constr in constraints
        if typeof(constr) <: Function
            push!(cbfuncs, constr)
            push!(constrtypes, :PB)
        else
            push!(cbfuncs, constr[2])
            push!(constrtypes, constr[1])
        end
    end
    Evaluator(cbfuncs, constrtypes, ndims, mdims, 
              printcxx = printcxx, vectorin = vectorin)
end
function Evaluator(cbfuncs, constrainttypes, ndims, mdims; 
                   printcxx = false, vectorin = true, resbyref = false)
    CXXID.ev += 1
    cppclassname = "NOMAD_Evaluator_$(CXXID.ev)"
    cbfunc = createcbfunc(CXXID.ev, cbfuncs, ndims, vectorin, resbyref)
    argstring = join(["x[$i].value()" for i in 0:ndims-1], ",")
    cxxstring = """
        class $cppclassname : public NOMAD::Evaluator {
            public:
                $cppclassname  ( const NOMAD::Parameters & p ) : NOMAD::Evaluator ( p ) {}
        ~$cppclassname ( void ) {}
        bool eval_x ( NOMAD::Eval_Point   & x          ,
                      const NOMAD::Double & h_max      ,
                      bool                & count_eval  ) const
            {
            double *res = new double [$mdims];
            \$($(cbfunc))($argstring, res);
            for (int i = 0; i < $mdims; i++) x.set_bb_output  ( i, res[i]); 
            count_eval = true; 
            return true;
            }        
        };"""
    if printcxx; println(cxxstring); end
    eval(Cxx.process_cxx_string(cxxstring))
    Evaluator(ndims, mdims, cppclassname, cbfunc, constrainttypes)
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
    bboutdefstring = "bbot[0] = NOMAD::OBJ;"
    for (i, typ) in enumerate(ev.constrainttypes)
        bboutdefstring *= "\n        bbot[$i] = NOMAD::$typ;"
    end
    cxxstring = """int $funcname(double* bf_x, double* bf_f, 
                                 NOMAD::Stats *stats) {
        NOMAD::Display out ( std::cout );
        out.precision ( NOMAD::DISPLAY_PRECISION_STD );
        try {
            NOMAD::Parameters p (out); 
            p.set_DIMENSION ($(ev.ndims));
            vector<NOMAD::bb_output_type> bbot ($(ev.mdims));
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

function optimize(opt)
    bf_x = zeros(opt.evaluator.ndims)
    bf_f = [0.]
    stats = icxx"new NOMAD::Stats();"
    opt.cppfunc(bf_x, bf_f, stats)
    bf_x, bf_f[1], icxx"$stats->get_bb_eval();"
end
export optimize

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
    params::Array{Tuple{Symbol, Any}, 1}
end
NOMADSolver() = NOMADSolver([(:DISPLAY_DEGREE, 0)])

function SolverInterface.NonlinearModel(s::NOMADSolver)
    NOMADMathProgModel(nothing, nothing, :Min, Float64[], 
                       s.params, 0., :SUCCESS)
end

mutable struct NOMADMathProgModel <: SolverInterface.AbstractNonlinearModel
    ev
    opt
    sense::Symbol
    x::Array{Float64, 1}
    params::Array{Tuple{Symbol, Any}, 1}
    obj::Float64
    status::Symbol
end

replaceinf(b) = map(x -> isinf(x) ? "$(sign(x) == 1 ? "" : "-")INFINITY" : x, b)

function getboundfunc(b, i, sign, d, numVar)
    constrval = Array{Float64}(numVar)
    function g(x)
        SolverInterface.eval_g(d, constrval, x)
        sign * (b - constrval[i])
    end
end

function MathProgBase.loadproblem!(m::NOMADMathProgModel, numVar, numConstr, 
                                   l, u, lb, ub, sense,
                                   d::SolverInterface.AbstractNLPEvaluator)
    SolverInterface.initialize(d, Symbol[])
    mdims = 1
    for i in 1:numConstr
        if !isinf(lb[i]); mdims += 1; end
        if !isinf(ub[i]); mdims += 1; end
    end
    function cbfunc(x, res)
        obj = sense == :Max ? -SolverInterface.eval_f(d, x) :
                                    SolverInterface.eval_f(d, x)
        unsafe_store!(res, obj, 1)
        if numConstr > 0
            constrval = Array{Float64}(numConstr)
            SolverInterface.eval_g(d, constrval, x)
            k = 1
            for i in 1:numConstr
                if !isinf(lb[i])
                    unsafe_store!(res, lb[i] - constrval[i], k += 1)
                end
                if !isinf(ub[i])
                    unsafe_store!(res, constrval[i] - ub[i], k += 1)
                end
            end
        end
    end
    m.ev = Evaluator([cbfunc], fill(:EB, mdims - 1), numVar, mdims, resbyref = true)
    m.params = [m.params; (:UPPER_BOUND, replaceinf(u)); (:LOWER_BOUND, replaceinf(l))]
end
function SolverInterface.setwarmstart!(m::NOMADMathProgModel,x)
    push!(m.params, (:X0, copy(float(x))))
end

function SolverInterface.optimize!(m::NOMADMathProgModel)
    isa(m.ev, Evaluator) || error("Must load problem before solving")
    m.opt = Opt(m.ev; m.params...)
    m.x, m.obj, nbbcalls = optimize(m.opt)
end

function SolverInterface.status(m::NOMADMathProgModel)
    if m.status == :SUCCESS || m.status == :FTOL_REACHED || m.status == :XTOL_REACHED
        return :Optimal
    elseif m.status == :ROUNDOFF_LIMITED
        return :Suboptimal
    elseif m.status in (:STOPVAL_REACHED,:MAXEVAL_REACHED,:MAXTIME_REACHED)
        return :UserLimit
    else
        error("Unknown status $(m.status)")
    end
end

SolverInterface.getsolution(m::NOMADMathProgModel) = m.x
SolverInterface.getobjval(m::NOMADMathProgModel) = m.sense == :Max ? -m.obj : m.obj

end # module
