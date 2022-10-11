# ovCancer-signature - Wilcoxon filter

#  Copyright (C) 2022 Fernando Palluzzi
#  e-mail: <fernando.palluzzi@gmail.com>
#  Bioinformatics facility 
#  Gemelli Science and Technological Park (GSTeP)
#  Fondazione Policlinico Universitario Agostino Gemelli IRCCS,
#  Largo Agostino Gemelli 8, 00168 Roma, Italy

#  ovCancer-signature is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  ovCancer-signature is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.


W <- data.frame(symbol = vector(), W = vector(), estimate = vector(), ci95 = vector(), pvalue = vector())
for (j in 2:43) {
	w <- wilcox.test(data[1:25, j], data[26:44, j], conf.int = TRUE)
	W <- rbind(W, data.frame(symbol = colnames(data)[j],
	                         W = w$statistic,
	                         estimate = w$estimate,
	                         ci95 = paste0(w$conf.int, collapse = ", "),
	                         pvalue = w$p.value))
}

top <- W$symbol[W$pvalue < 0.05]
top
