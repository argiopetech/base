select "Multi-Pop" as program, runId, time, count(iterId)
  from run
    join multi_pop
  on run.id = runId
  group by runId
union
select "Single Pop", runId, time, count(iterId)
  from run
    join single_pop
  on run.id = runId
  group by runId
order by runId;
