using SimpleDelayed

task_id = parse(Int, ARGS[1])
println("Running task $task_id")
# TODO: map task_id to (system, LFMO model, delay option),
# run the experiment, write output/result_$(task_id).csv
