select h.start_longitude, h.start_latitude, h.bottom_depth, h.stationid, c.* 
from racebase.catch c, racebase.haul h 
where h.stationid in ('61-15', '61-02', '11-34', '11-29', '11-41', '11-18', 'A-02', 'AZ0504') 
and c.species_code in (401, 402, 403, 411, 421, 436, 441, 446, 456, 461, 473, 474, 476, 478, 481, 484, 486)
and c.hauljoin = h.hauljoin 
order by h.stationid, c.number_fish, c.cruise desc


select h.start_longitude, h.start_latitude, h.bottom_depth, h.stationid, c.* 
from racebase.catch c, racebase.haul h 
where h.stationid in ('61-15', '61-02', '11-34', '11-29', '11-41', '11-18', 'A-02', 'AZ0504') 
and c.species_code in (401, 402, 403, 411, 421, 436, 441, 446, 456, 461, 473, 474, 476, 478, 481, 484, 486)
and c.hauljoin = h.hauljoin 
order by c.weight desc



select h.*
from racebase.haul h 
where h.stationid in ('61-15', '61-02', '11-34', '11-29', '11-41', '11-18', 'A-02', 'AZ0504')  
and h.performance >= 0