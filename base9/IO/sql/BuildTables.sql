R"=====(
create table if not exists run
( id   integer not null
, time datetime default CURRENT_TIMESTAMP
, primary key (id));

create table if not exists iteration
( runId integer not null
, id    integer not null
, primary key (runId, id)
, foreign key (runId) references run);

create table if not exists single_pop
( runId          integer not null
, iterId         integer not null
, logAge         double
, y              double
, feh            double  not null
, modulus        double  not null
, absorption     double  not null
, carbonFraction double
, ifmrIntercept  double
, ifmrSlope      double
, ifmrQuadCoef   double
, logPost        double  not null
, stage          integer not null
, primary key (runId, iterId)
, foreign key (runId, iterId) references iteration);

create table if not exists multi_pop
( runId          integer not null
, iterId         integer not null
, logAge         double
, feh            double  not null
, modulus        double  not null
, absorption     double  not null
, carbonFraction double
, ya             double  not null
, yb             double  not null
, lambda         double  not null
, logPost        double  not null
, stage          integer not null
, primary key (runId, iterId)
, foreign key (runId, iterId) references iteration);

create table if not exists field_star_likelihood
( runId integer not null
, likelihood double not null
, primary key (runId)
, foreign key (runId) references run(id));

create table if not exists star
( runId integer not null
, starId text not null
, primaryMass double not null
, secondaryMassRatio double not null
, stage integer not null
, clusterMembershipPrior double not null
, useDuringBurnin boolean not null
, primary key (runId, starId)
, foreign key (runId) references run);

create table if not exists photometry
( runId integer not null
, starId text not null
, filter text not null
, magnitude double not null
, stdDev double not null
, foreign key (runId, starId) references star);

create table if not exists star_posterior
( runId integer not null
, iterId integer not null
, starId text not null
, clust_a_posterior double not null
, clust_b_posterior double not null
, primary key (runId, iterId, starId)
, foreign key (runId, iterId) references iteration (runId, id)
, foreign key (runId, starId) references star (runId, starId));

create table if not exists sample_mass
( runId integer not null
, referencedRunId integer not null
, iterId integer not null
, starId text not null
, primaryMass double not null
, massRatio double not null
, clusterMembership double not null
, primary key (runId, referencedRunId, iterId, starId)
, foreign key (referencedRunId, iterId) references iteration (runId, id)
, foreign key (referencedRunId, starId) references star (runId, starId)
, foreign key (referencedRunId, iterId) references single_pop (runId, iterId));

create table if not exists sample_wd_mass
( runId integer not null
, referencedRunId integer not null
, iterId integer not null
, starId text not null
, mass double not null
, clusterMembership double not null
, precursorLogAge double not null
, coolingAge double -- Can sometimes be NaN, which registers as null
, logTeff double not null
, logG double not null
, primary key (runId, referencedRunId, iterId, starId)
, foreign key (referencedRunId, iterId) references iteration (runId, id)
, foreign key (referencedRunId, starId) references star (runId, starId)
, foreign key (referencedRunId, iterId) references single_pop (runId, iterId));
)====="
