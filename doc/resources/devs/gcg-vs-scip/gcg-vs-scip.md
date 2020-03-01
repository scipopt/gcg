# Differences between GCG and SCIP {#gcg-vs-scip}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

@todo add fancy table here

<table>
  <tr>
    <th>GCG</th>
    <th>SCIP</th>
  </tr>
  <tr>
    <td>
    <div class="fragment">
      <div class="line">SCIP Status        : problem is solved [optimal solution found]</div>
      <div class="line">Solving Time (sec) : 28.48</div>
      <div class="line">Solving Nodes      : 1</div>
      <div class="line">Primal Bound       : -4.10000000000000e+01 (5 solutions)</div>
      <div class="line">Dual Bound         : -4.10000000000000e+01</div>
      <div class="line">Gap                : 0.00 %</div>
    </div>
    </td>
    <td>
    <div class="fragment">
      <div class="line">SCIP Status        : problem is solved [optimal solution found]</div>
      <div class="line">Solving Time (sec) : 337.72</div>
      <div class="line">Solving Nodes      : 1088063 (total of 1088983 nodes in 2 runs)</div>
      <div class="line">Primal Bound       : -4.10000000000000e+01 (389 solutions)</div>
      <div class="line">Dual Bound         : -4.10000000000000e+01</div>
      <div class="line">Gap                : 0.00 %</div>
    </div>
    </td>
  </tr>
</table>

GCG allows much easier implementation of anything that is somehow related to branch-and-price, since it has its own universal toolbox of methods.

It may facilitate the understanding to know about the interaction between GCG and SCIP:

@subpage mirroring "Interaction during Branching"
@subpage preprocessing
