# Known Bugs {#known-bugs}

<table>
  <tr>
    <th>Bug</th>
    <th>Reason</th>
    <th>Version</th>
  </tr>
  <tr>
    <td>Pressing CTRL-C does not always terminate GCG</td>
    <td>GCG creates SCIP instances for subproblems. If CTRL-C is pressed, one of these subinstances is stopped, but not the calling function of GCG.<br></td>
    <td>3.0.0, current master</td>
  </tr>
</table>
