using Printf, JSON
include("pipeline2b.jl")
using .Pipeline2b
using .Pipeline2b.MomentsVerified
nodes,wts = MomentsVerified.gl_reference(24)
println("GL reference certified (n_p=24)."); flush(stdout)
@printf("\n===== Omega2 Stage-A N=64 (lean assembler) =====\n"); flush(stdout)
rec = Pipeline2b.run_N(64; npanel=96, nodes=nodes, wts=wts)
open("stageA_dom2_N64.json","w") do f; JSON.print(f,[rec],2); end
@printf("  [checkpoint] N=64 bracket=[%.10f,%.10f] L2=[%.12f,%.12f] sep=%s\n",
  rec["bracket_lo"],rec["bracket_hi"],rec["L2"][1],rec["L2"][2],rec["sep_L2_gt_U1"]); flush(stdout)
println("STAGEA_N64_DONE")
