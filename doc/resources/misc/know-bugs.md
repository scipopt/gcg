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
  <tr>
    <td>As described in the SCIP documentation, one should be able to run a test with
    `SETTINGS="set1,set2"` and SCIP should run the test on both settings. This does not work.</td>
    <td>Unknown</td>
    <td>3.0.0, current master</td>
  </tr>
</table>
