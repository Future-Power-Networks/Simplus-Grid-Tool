# Debug Logs

* 10-03-2020, frame transformation

    Investigate the effect of frame transformatoins. The result is that in-file transformations are not as stable as out-file (in simulink blocks) and the reasons are mysterious and not revealed yet (see commit test_transform, d69d5a5). If the input transformation is taken out, it becomes stable (see commit half_half, c6cc5b6).
