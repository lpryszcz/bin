#!/usr/bin/bash
# pip install spotdl
outdir=~/Music/spotify/spotdl
echo "Saving files to $outdir ..."
mkdir -p $outdir
# download daily 1-6 & release radar
cd $outdir
args="--threads 6 --overwrite skip"
spotdl $ags --m3u ../`date +%F`.DailyMix1.m3u https://open.spotify.com/playlist/37i9dQZF1E3783ZEsHDeht
sed -i 's,^spotdl/,,;s,^,spotdl/,' ../`date +%F`.*.m3u
spotdl $ags --m3u ../`date +%F`.DailyMix2.m3u https://open.spotify.com/playlist/37i9dQZF1E38RsxzxUTZRk
spotdl $ags --m3u ../`date +%F`.DailyMix3.m3u https://open.spotify.com/playlist/37i9dQZF1E38cTYLyGroCC
spotdl $ags --m3u ../`date +%F`.DailyMix4.m3u https://open.spotify.com/playlist/37i9dQZF1E36VpRu7l05xK
spotdl $ags --m3u ../`date +%F`.DailyMix5.m3u https://open.spotify.com/playlist/37i9dQZF1E390FOdL3Y9Gq
spotdl $ags --m3u ../`date +%F`.DailyMix6.m3u https://open.spotify.com/playlist/37i9dQZF1E38dIUgR6Sw1c
spotdl $ags --m3u ../`date +%F`.ReleaseRadar.m3u https://open.spotify.com/playlist/37i9dQZEVXbv5k5KXCSJOI
spotdl $ags --m3u ../`date +%F`.DiscoverWeekly.m3u https://open.spotify.com/playlist/37i9dQZEVXcJgVGv6YnNez

# add spotdl to m3u (strip first if added before)
sed -i 's,^spotdl/,,;s,^,spotdl/,' ../`date +%F`.*.m3u


# other playlists
#spotdl $ags --m3u ../KitchenSwagger.m3u https://open.spotify.com/playlist/37i9dQZF1DX2FsCLsHeMrM
#spotdl $ags --m3u ../PiscoLoungeBar.m3u https://open.spotify.com/playlist/6Qm5HyjYBV9NCBqHWWLt04
#spotdl $ags --m3u ../NastyMondaysApolo.m3u https://open.spotify.com/playlist/0Rc1bR5YLZ8fKXAZFiwkPd
#spotdl $ags --m3u ../FIB2020.m3u https://open.spotify.com/playlist/4QWtWxXa8Cc4RdqQTP7BV2
#spotdl $ags --m3u ../VanLife.m3u https://open.spotify.com/playlist/37i9dQZF1DX2ogDiL6nZJr
#spotdl $ags --m3u ../TotallyStressFree.m3u https://open.spotify.com/playlist/37i9dQZF1DWT7XSlwvR1ar
#spotdl $ags --m3u ../.m3u 
#spotdl $ags --m3u ../.m3u
# ulimit -n 40000 # you may need to edit /etc/security/limits.conf
#spotdl $ags --m3u ../LikedSongs.m3u https://open.spotify.com/playlist/1jy5hPElIVzLBxryskljkc #saved --user-auth --client-id 69229da242084239b102dc3d113e1337 --client-secret dc69f1749c124de39e1782ad2b8b9034
#sed -i 's,^spotdl/,,;s,^,spotdl/,' ../*.m3u
