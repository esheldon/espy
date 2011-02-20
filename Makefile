#  This file is part of Erin Sheldon's personal IDL routines.
#
#  Copyright (C) 2005  Erin Sheldon, NYU.  erin.sheldon at gmail.com
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


PREFIX := ""

copy :
ifeq ($(PREFIX),"")
	@echo PREFIX is not set, aborting
else

	@echo "rsyncing to $(PREFIX)"
	@if [ ! -e $(PREFIX) ]; then \
		echo Creating directory $(PREFIX); \
		mkdir -p $(PREFIX); \
	fi
	rsync -av \
		--exclude "*svn*" \
		--exclude "*git*" \
		--exclude "*swp" \
		--exclude "*~" \
		--exclude "*pyc" ./ $(PREFIX)/

endif


