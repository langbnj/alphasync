CREATE TABLE `alphacon` (
  `acc` char(13) NOT NULL,
  `site1` mediumint NOT NULL,
  `site2` mediumint NOT NULL,
  `aa1` char(1) DEFAULT NULL,
  `aa2` char(1) DEFAULT NULL,
  `atom1` char(3) NOT NULL,
  `atom2` char(3) NOT NULL,
  `type` enum('AromaticContacts','CarbonylContacts','CovalentContacts','HBondContacts','HydrophobicContacts','IonicContacts','MetalContacts','PolarHBondContacts','VanDerWaalsContacts','WeakHBondContacts','WeakPolarHBondContacts') NOT NULL,
  `afdb` char(1) NOT NULL,
  `dist` float DEFAULT NULL,
  `pae` float DEFAULT NULL,
  PRIMARY KEY (`acc`,`site1`,`site2`,`atom1`,`atom2`,`type`,`afdb`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='AlphaSync contacts';

CREATE TABLE `alphafrag` (
  `acc` varchar(15) NOT NULL,
  `name` varchar(16) DEFAULT NULL,
  `species` varchar(5) DEFAULT NULL,
  `tax` mediumint DEFAULT NULL,
  `frag` smallint NOT NULL,
  `fragstart` mediumint DEFAULT NULL,
  `fragstop` mediumint DEFAULT NULL,
  `source` varchar(30) NOT NULL,
  `afdb` tinyint NOT NULL,
  `seq` varchar(36000) DEFAULT NULL,
  PRIMARY KEY (`acc`,`frag`,`afdb`,`source`),
  KEY `Name` (`name`),
  KEY `Species` (`species`),
  KEY `Frag` (`frag`),
  KEY `Fragstart` (`fragstart`),
  KEY `Fragstop` (`fragstop`),
  KEY `Source` (`source`),
  KEY `Alphafold_db` (`afdb`),
  KEY `Seq` (`seq`(50)),
  KEY `Tax` (`tax`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='AlphaFold fragment sequences';

CREATE TABLE `alphamap` (
  `type` char(12) NOT NULL,
  `version` char(7) NOT NULL,
  `value` char(23) NOT NULL,
  `best` tinyint NOT NULL,
  `afdb` tinyint NOT NULL,
  `map` char(13) NOT NULL,
  `species` char(32) DEFAULT NULL,
  `tax` mediumint DEFAULT NULL,
  `avg_plddt` float DEFAULT NULL,
  PRIMARY KEY (`type`,`version`,`value`,`best`,`afdb`,`map`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='AlphaSync mapping from latest UniProt to AlphaFold DB accessions';

CREATE TABLE `alphasa` (
  `acc` char(13) NOT NULL,
  `site` mediumint NOT NULL,
  `afdb` char(1) NOT NULL,
  `aa` char(1) DEFAULT NULL,
  `plddt` float DEFAULT NULL,
  `plddt10` float DEFAULT NULL,
  `asa` float DEFAULT NULL,
  `asa10` float DEFAULT NULL,
  `relasa` float DEFAULT NULL,
  `relasa10` float DEFAULT NULL,
  `dis` char(1) DEFAULT NULL,
  `dis10` char(1) DEFAULT NULL,
  `surf` char(1) DEFAULT NULL,
  `surf10` char(1) DEFAULT NULL,
  `sec` char(1) DEFAULT NULL,
  `iso` char(1) DEFAULT NULL,
  `phi` float DEFAULT NULL,
  `psi` float DEFAULT NULL,
  `omega` float DEFAULT NULL,
  `chi1` float DEFAULT NULL,
  `chi2` float DEFAULT NULL,
  `chi3` float DEFAULT NULL,
  `chi4` float DEFAULT NULL,
  `chi5` float DEFAULT NULL,
  `tau` float DEFAULT NULL,
  PRIMARY KEY (`acc`,`site`,`afdb`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='AlphaSync accessible surface areas';

CREATE TABLE `alphaseq` (
  `acc` char(13) NOT NULL,
  `name` char(16) DEFAULT NULL,
  `species` char(5) DEFAULT NULL,
  `tax` mediumint DEFAULT NULL,
  `frags` smallint DEFAULT NULL,
  `afdb` tinyint NOT NULL,
  `seq` varchar(36000) DEFAULT NULL,
  PRIMARY KEY (`acc`,`afdb`),
  KEY `Name` (`name`),
  KEY `Species` (`species`),
  KEY `Seq` (`seq`(50)),
  KEY `Tax` (`tax`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='AlphaFold protein sequences';

CREATE TABLE `alphastats` (
  `stat` varchar(50) NOT NULL,
  `value` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`stat`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='AlphaSync precalculated statistics';

CREATE TABLE `alphauniprot` (
  `acc` varchar(13) NOT NULL,
  `canon` varchar(10) DEFAULT NULL,
  `name` varchar(16) DEFAULT NULL,
  `species` varchar(20) DEFAULT NULL,
  `tax` mediumint DEFAULT NULL,
  `species_common` varchar(250) DEFAULT NULL,
  `species_latin` varchar(250) DEFAULT NULL,
  `reviewed` tinyint DEFAULT NULL,
  `refproteome` tinyint DEFAULT NULL,
  `refprotcanon` tinyint DEFAULT NULL,
  `fullname` varchar(1000) DEFAULT NULL,
  `symbols` varchar(1000) DEFAULT NULL,
  `synonyms` varchar(1000) DEFAULT NULL,
  `func` varchar(15000) DEFAULT NULL,
  `seqlen` mediumint DEFAULT NULL,
  `seq` varchar(46000) DEFAULT NULL,
  PRIMARY KEY (`acc`),
  KEY `Canon` (`canon`),
  KEY `Name` (`name`),
  KEY `Species` (`species`),
  KEY `Seq` (`seq`(50)),
  KEY `Tax` (`tax`),
  KEY `Symbols` (`symbols`),
  KEY `Synonyms` (`synonyms`),
  KEY `Fullname` (`fullname`),
  KEY `Reviewed` (`reviewed`),
  KEY `Refproteome` (`refproteome`),
  KEY `Refprotcanon` (`refprotcanon`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='UniProt annotation via API';

CREATE TABLE `alphauniprot_species` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `tax` mediumint NOT NULL,
  `species` varchar(20) DEFAULT NULL,
  `species_latin` varchar(250) DEFAULT NULL,
  `species_latin_short` varchar(250) DEFAULT NULL,
  `species_common` varchar(250) DEFAULT NULL,
  `fullspecies` varchar(250) DEFAULT NULL,
  `compara_species` tinyint DEFAULT NULL,
  `complete` tinyint DEFAULT NULL,
  `complete_iso` tinyint DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Tax` (`tax`),
  KEY `Species` (`species`),
  KEY `Species_latin` (`species_latin`),
  KEY `Species_latin_short` (`species_latin_short`),
  KEY `Species_common` (`species_common`),
  KEY `Fullspecies` (`fullspecies`),
  KEY `compara_species` (`compara_species`),
  KEY `Complete` (`complete`)
  KEY `Complete_iso` (`complete_iso`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='UniProt species summary table';

CREATE TABLE `alphauniprot_symbols` (
  `acc` char(13) NOT NULL,
  `species` varchar(20) DEFAULT NULL,
  `tax` mediumint NOT NULL,
  `type` enum('symbol','synonym') NOT NULL,
  `alias` varchar(500) NOT NULL,
  PRIMARY KEY (`alias`,`type`,`tax`,`acc`),
  KEY `Species` (`species`),
  KEY `Type` (`type`),
  KEY `Tax` (`tax`),
  KEY `Acc` (`acc`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='UniProt normalised gene symbols and synonyms';

CREATE TABLE `alphasync_api_log`.`rate_limits` (
  `hashed_ip` char(64) NOT NULL,
  `requests` int DEFAULT NULL,
  `timestamp` int DEFAULT NULL,
  PRIMARY KEY (`hashed_ip`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
