create table if not exists public.model_param(
    "parID" int not null,
    "circuit_n" text not null ,
    "variant" text not null,
    "n_samples" int not null,
    primary key ("parID", "circuit_n", "variant", "n_samples")
);

-- create table if not exists analytical_result(
--     "parID" int not null
-- )



