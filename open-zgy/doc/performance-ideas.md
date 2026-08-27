## Notes about OpenZGY performance

The speed of OpenZGY and of the old ZGY-Public library mostly depends
on the cost of the raw I/O. The performance of the new and the old
code is expected to be similar. This is supported by fair benchmarks.
So far, there is no smoking gun that shows otherwise.

OpenZGY is better than the old code when it comes to compression and
decompression. And, depending on access patterns, for cloud I/O where
the old code supported that at all. Probably (hopefully) OpenZGY has
many other places where it is faster and/or better than the old code.
This document, however, focuses on the opposite situation.

Some benchmarks done in lab conditions have indicated that OpenZGY is
slower. The ones reported to date have done very little I/O. E.g.,
because the files were small enough to fit in memory. And possibly
being pre-loaded due to the test being run several times.

Some tests completely do away with I/O. E.g., reading a Petrel
calculator cube instead of real data. Such tests are expected to run
very fast. Several of them are shown to run slower, but still quite
fast, in OpenZGY than they do in ZGY-Public. The test results will
often not be relevant to ordinary use of the application.

The above doesn't mean that the OpenZGY team is sure there are no
noticeable performance problems. In some cases there might be a
trade-off between performance and accuracy. Other cases might be due to
unfortunate design. Or even straight out bugs.

## Details, for OpenZGY developers only.

Here is a list of tasks that might be considered that could
potentially improve performance. With stress on potentially. Please
don't judge the OpenZGY code by the size of this list. And don't waste
time implementing any of the suggestions without at least trying to
measure whether it would help. Most won't.

### Cloud I/O

* Implement a large cache, or implement "rema-1000" caching of just a
  couple of bricks, to work around applications reading just one brick
  at a time. This would help Petrel. But it would help a lot more to
  fix this inside Petrel. Which does in fact have its own cache. It
  just isn't using it optimally.

* Parallel cloud I/O and CPU tasks. Even when CPU load is much
  smaller that I/O, this might in theory have an effect.

* Change IOContext defaults to what is most likely the fastest.
  Instead of aiming for conservative resource usage. But, it is much
  better to have the application tailor this to its own access
  pattern.

* Allow an application to provide a different IOContext for finalize()
  than it does for write(). This allows for more tweaking by the
  application. It doesn't change anything by itself.

### Local I/O

* On windows, the native ReadFile / WriteFile might be faster than
  Posix I/O. Could even try asynchronous I/O. The old accessor does
  this. Initial measurements don't show much difference, but it might
  depend on usage pattern.

* Pre-allocate space on local files to allow writes to be done in
  parallel. Needs help from application to know which bricks will not
  end up all constant. Probably not much gain, and is only applicable
  to local files, and would place a large burden on the application.
  Petrel will not benefit without significant changes.

### Both local and cloud I/O

* Incremental low-res computation, possibly only for LOD 1. High cost.
  The old accessor does this. Measurements show a significant
  improvement but only in some scenarios. Typically when files are
  very large. In particular, large files output to cloud.

* Parallelize more of lod compute, with 4 threads (only). Note that
  increased memory consumption might kill any benefits.

* Short cut to generate LOD1 data immediately if the application
  chooses to write either 8 bricks or 4 brick columns at a time.
  Petrel will not benefit.

### CPU-bound

* A new brick API and a new implementation of genlod making use of
  that should speed up finalize. There will be less buffer copies and
  a smaller RAM working set. Down to 9 MB, which fits comfortably into
  L3$ which us typically 16-64 MB. The changes inside finalize will be
  expensive. Also, need to keep the old algorithm for use in LOD0 and all
  cloud LODs. Also, when multi-threading the L3$ is shared between
  cores. So it won't be big enough anyway.

* A new optional brick API might also help CPU usage of some
  application read/write. Applications opting in would need to change.
  Applications that cannot change but currently read/write just one
  brick at a time might see a slight performance boost. Because the
  brick API might be able to utilize more shortcuts.
  So, Petrel might see a small performance boost.

* int->float, reshape, and range() done in the same loop to reduce
  OpenMP overhead. Need to actually implement this before I can
  measure the effect. But, probably not much.

* Carefully measure that no OpenMP loops add so much overhead that
  they end up running slower. Already done, but I could try even
  harder. The penalty for having too many OpenMP loops can be a lot
  worse that having too few. Specifically, check the OpenMP loop
  calling RoundAndClip for 256*1024 samples.

* CopySubset to recognize more types of shortcut. Lower cost and risk
  but I doubt it will help much. The cost might lie in warming the CPU
  cache. Likely N/A if implementing and using a brick API.

* Profile and optimize LodWeightedAverage; the main CPU hog apart from
  compression. OR, replace it with a cheaper, lower quality algorithm.
  Not implemented in the old accessor. On my machine, total CPU time
  went down 10% in the writeonlly/smallfile/singlethread case. Much
  less when multi-threaded. Caveat, those MT measurements are a bit
  flaky. Still, it looks like this is not a significant help.
  If it does in fact make a difference, then Petrel would benefit.

* Special shortcut when reading or writing a single brick. Which is
  discouraged, but Petrel is going to do this anyway.
  Limit the shortcut to no sdms, consolidate-writes, compression,
  byte swapping, conversion to storage, read/modify/write.
  This is a simpler version of the proposed brick api.
  So, Petrel might see a small performance boost. Others probably not.

* As the single-brick shortcut, but also require that the application
  cannot mix the brick api and the general api. implement a different
  threading model where the application and not the OpenZGY library is
  doing the threading. Petrel might see a significant speedup in
  write. But not in finalize.

* Allow the application to tell how much memory OpenZGY is allowed to
  use for finalize. Instead of using a conservative estimate. Fairly
  low cost and risk, but also fairly low expectations.

* Lazy evaluation of apparent block size and brick-end. Leads to less
  consistency checking. Petrel benefits in that opening a file might
  be faster. But overall time would remain the same.

* Caching some derived metadata. Very unlikely to make a difference.

* Have genlod do its own memory management and re-use read buffers.
  Unlikely to help much.

* Skip initializing buffers to zero. The effect is questionable and
  the benefits are probably marginal. Note that the sampling profiler
  complains that the operation is expensive. But the actual cost seems
  to be moving data in and out of the (probably L3) cache. By skipping
  the initialize the cost is just moved to the next operation. Note
  that in a future brick-API the change is likely safe. So, keep it in
  this list.

* Several small functions might be rewritten to use sse2 intrinsics.
  See [doc/sse2.md](sse2.md) for details.

## Drawbacks

Several of the suggestions relate to some kind of trade-offs. This
means that several issues in the old ZGY-Public might come back.
Many of the issues below are caused by incremental lod generation.

* The contents of the low resolution data might show visible brick
  artifacts due to the low-pass filter in LOD1 not being fed entire
  traces. (speed / quality trade-off).

* There might be additional brick artifacts in LOD2 and higher due to
  the weighted average algorithm. Especially between the first brick
  written and the first brick written with "typical" sample values for
  this data set. Might mitigate this by using a different algorithm
  (choosing a different kind of lower quality), or only do LOD1
  incrementally (speed / quality trade-off).

* The low resolution data will be significantly noisier when
  compression is enabled, and gets noisier for higher lod levels. Past
  a certain point the low resolution data might be completely useless.
  Can be partly mitigated by not compressing higher levels of lods
  (space / quality trade-off).

* The histogram range will typically end up wider than needed, so the
  effective number of bins might be even lower than 256. Can be partly
  mitigated by using a wider histogram internally. As a side benefit,
  the need to track statistics for int8 and int16 would go away. Note
  that even with the mitigation, there would still be odd cases where
  the histogram ends up useless.

* The statistical range might be wrong if a spike was overwritten. It
  will be wrong if that spike was the one and only sample responsible
  for either the min or the max statistical range.

* Both the histogram range, the statistical range, and the low
  resolution data might be slightly non-deterministic. Depending on
  the order the data is written and/or updated.

* If implementing the short cut to generate LOD1 data immediately, and
  compression is enabled, the problems with non deterministic lod
  levels and with low resolution brick artifacts become significantly
  worse.

* Files that were written to the cloud with incremental lod generation
  might be slightly slower to read later. This is because of fewer
  contiguous bricks when low-res bricks are interleaved. (write speed /
  read speed trade-off). Partial mitigation to only compute low-res when
  4 entire brick columns (instead of 8 neighboring bricks) have been
  written.

* Generating lod incrementally when the output is on the cloud will in
  some cases be expensive due to reading too small blocks. Mitigation
  as above; compute full brick-columns.

## Benefits other than possibly better performance

* When not generating lod incrementally, i.e. with today's
  implementation, it is a bit tricky to show a progress report to the
  user. If the application cannot use the progress callback in
  finalize() then an existing progress bar might spend approximately
  the same time going from 0%-99% as it does going from 99% to 100%.
  So it appears to be stuck at 99%. This is a problem when realizing
  to ZGY in Petrel. I don't know how easy it is to fix inside Petrel.
  Fixing it in OpenZGY can only be done with incremental lod
  generation. Which will trigger almost all the drawbacks listed
  above.
