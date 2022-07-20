#!/usr/bin/bash
outdir=~/Music/spotify/`date +%F`
echo "Saving files to $outdir ..."
mkdir -p $outdir
# download daily 1-6, discover weekly & release radar
cd $outdir && spotify_dl -s -w -k -l \
	   https://open.spotify.com/playlist/37i9dQZF1E3783ZEsHDeht \
	   https://open.spotify.com/playlist/37i9dQZF1E38RsxzxUTZRk \
	   https://open.spotify.com/playlist/37i9dQZF1E38cTYLyGroCC \
	   https://open.spotify.com/playlist/37i9dQZF1E36VpRu7l05xK \
	   https://open.spotify.com/playlist/37i9dQZF1E390FOdL3Y9Gq \
	   https://open.spotify.com/playlist/37i9dQZF1E38dIUgR6Sw1c \
	   https://open.spotify.com/playlist/37i9dQZEVXcJgVGv6YnNez \
	   https://open.spotify.com/playlist/37i9dQZEVXbv5k5KXCSJOI


