import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.nddata import NDData
from astropy.visualization import ZScaleInterval
from astroquery.gaia import Gaia
from lmfit import Model, Parameters
from photutils.isophote import Ellipse, EllipseGeometry
from photutils.psf import extract_stars
from photutils.centroids import centroid_2dg
from photutils.detection import IRAFStarFinder
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


def detect_band_name(header, fits_path):
    """
    Dynamically detects whether the UVIT image belongs to FUV (starts with F)
    or NUV (starts with N) checking header keywords, file name, and parent folder path.
    """
    fname = os.path.basename(fits_path).upper()
    parent_dir = os.path.basename(os.path.dirname(fits_path)).upper()
    filter_val = str(header.get('FILTER', '')).upper().strip()
    band_val = str(header.get('BAND', '')).upper().strip()
    candidates = [filter_val, band_val, fname, parent_dir]
    for item in candidates:
        if not item:
            continue
        tokens = item.replace('_', ' ').replace('-', ' ').replace('.', ' ').replace('=', ' ').split()
        for token in tokens:
            if token.startswith('F'):
                return token  # FUV band (e.g., F148W, F154W, F169M, F172M, FUV)
            elif token.startswith('N'):
                return token  # NUV band (e.g., N242W, N219M, N245M, N263M, N279N, NUV)
    return 'F148W'  # Fallback default


def process_uvit_image(fits_path, output_dir="output", cutout_size=81, max_candidates=15):
    """
    Automated PSF extraction pipeline specifically built for UVIT (FUV / NUV) images.
    """
    os.makedirs(output_dir, exist_ok=True)
    pxs = 0.417  # UVIT pixel scale (arcsec/pixel)

    print(f"\n[+] Loading UVIT FITS File: {fits_path}")
    with fits.open(fits_path) as hdul:
        img_hdr = hdul[0].header
        img_data = hdul[0].data.astype(float)
    band_name = detect_band_name(img_hdr, fits_path)
    print(f"[+] Detected Filter Band: UVIT-{band_name}")
    wcs = WCS(img_hdr)
    height, width = img_data.shape
    xc, yc = width // 2, height // 2

    # Query Gaia DR3 Catalog
    print("[+] Querying Gaia DR3 Catalog...")
    center_sky = pixel_to_skycoord(xc, yc, wcs=wcs)
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    Gaia.ROW_LIMIT = -1
    gaia_query = Gaia.cone_search_async(center_sky, radius=u.Quantity(14, u.arcmin)).get_results()
    gaia_coords = SkyCoord(gaia_query['ra'], gaia_query['dec'], unit='deg')
    gaia_query['xpix'], gaia_query['ypix'] = wcs.world_to_pixel(gaia_coords)
    xy_coords = np.column_stack((gaia_query['xpix'], gaia_query['ypix']))

    # Automatic background estimation if SKY keyword is missing
    sky_bkg = img_hdr.get('SKY', np.nanmedian(img_data))
    threshold = 20 * sky_bkg
    print(f"[+] Detecting candidate stars via IRAFStarFinder (Threshold: {threshold:.4f})...")
    iraf = IRAFStarFinder(xycoords=xy_coords,
                            fwhm=3.5,
                            threshold=threshold,
                            exclude_border=True,
                            roundness_range=(0.0, 2.0),
                            sharpness_range=(0.0, 2.0),
                            n_brightest=max_candidates)
    df = iraf(img_data)
    if df is None or len(df) == 0:
        print("[-] No valid star candidates found matching criteria.")
        return

    # Filter out stars in extreme outer edges (>2000 pixels from center)
    valid_indices = []
    for i, (x, y) in enumerate(zip(df['x_centroid'], df['y_centroid'])):
        dist_from_center = np.sqrt((x - xc)**2 + (y - yc)**2)
        if dist_from_center <= 2000:
            valid_indices.append(i)
    df = df[valid_indices]
    print(f"[+] Identified {len(df)} candidate stars.")

    # Show Candidate Overlay Map using Native WCSAxes
    _plot_field_overlay(img_data, df, xc, yc, img_hdr)

    # Extract Cutouts
    star_sky_coords = pixel_to_skycoord(df['x_centroid'], df['y_centroid'], wcs=wcs)
    table = Table()
    table['x'] = np.round(df['x_centroid'], 4)
    table['y'] = np.round(df['y_centroid'], 4)
    table['ra'] = np.round(star_sky_coords.ra, 4)
    table['dec'] = np.round(star_sky_coords.dec, 4)
    nddata = NDData(data=img_data)
    extracted_stars = extract_stars(nddata, table, size=cutout_size)

    # Process Candidates & Interactive Selection
    _fit_and_select_psfs(extracted_stars, table, sky_bkg, pxs, output_dir, img_hdr, band_name)


def _plot_field_overlay(data, df, xc, yc, header):
    """Displays overview map of detected UVIT stars using native Astropy WCSAxes."""
    wcs = WCS(header)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(data)
    im = ax.imshow(data, vmin=vmin, vmax=vmax, cmap='gray', origin='lower')
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Counts')
    # Candidate Markers
    ax.scatter(df['x_centroid'], df['y_centroid'], edgecolors='lime', facecolors='none', s=100, linewidths=2.0, label='Candidates')
    # Field Boundary (2000 px radius)
    circle = plt.Circle((xc, yc), 2000, color='red', fill=False, linewidth=1.5, label='FOV Limit')
    ax.add_patch(circle)
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.grid(color='white', alpha=0.3, linestyle='dotted')
    fig.suptitle("Detected UVIT Candidate Stars", fontsize=16, y=0.92)
    plt.show()


def _fit_and_select_psfs(stars, table, sky_bkg, pxs, output_dir, header, name):
    """Fits 1D Moffat+Exponential profiles and lets user inspect/select candidates interactively."""
    zscale = ZScaleInterval()
    accepted_candidates = []
    if 'SKY' not in header:
        header['SKY'] = sky_bkg
    for idx, star in enumerate(stars):
        idx = idx + 1
        data = star.data
        fname = f"Star_{idx}"
        print(f"\n==================== Candidate Star #{idx} ====================")
        
        # 1D Isophotal Fit
        c = centroid_2dg(data)
        geom = EllipseGeometry(x0=c[0], y0=c[1], sma=3, eps=0.01, pa=0.0)
        ellipse = Ellipse(data, geometry=geom)
        try:
            isolist = ellipse.fit_image(fix_eps=True, fix_center=True, minsma=0.01, maxsma=34.0, step=0.10, linear=False, minit=20000)
            psf_tab = Table(isolist.to_table())
            PSF = Table()
            PSF['sma_a'] = psf_tab['sma'] * pxs
            PSF['int_i'] = psf_tab['intens']
            PSF['int_i_e'] = psf_tab['intens_err']
            Im0 = PSF['int_i'][0]
            def exp_fn(r, Ie0, h0):
                return Ie0 * np.exp(-r / h0)
            def mof_fn(r, alpha, beta):
                return Im0 * (1 + (r / alpha)**2)**(-beta)
            def total_model(r, alpha, beta, Ie0, h0):
                return mof_fn(r, alpha, beta) + exp_fn(r, Ie0, h0)

            # Fit outer exponential tail
            psf_exp = PSF[PSF['sma_a'] > 6.0]
            par_exp = Parameters()
            par_exp.add_many(('Ie0', 5.0, True, None, None), ('h0', 2.0, True, None, None))
            res_exp = Model(exp_fn).fit(psf_exp['int_i'], par_exp, r=psf_exp['sma_a'], weights=1/psf_exp['int_i_e'])

            # Fit full profile
            par_full = Parameters()
            par_full.add_many(
                ('alpha', 1.5, True, None, None),
                ('beta', 1.5, True, None, None),
                ('Ie0', res_exp.params['Ie0'].value, False, None, None),
                ('h0', res_exp.params['h0'].value, False, None, None)
            )
            res_full = Model(total_model).fit(PSF['int_i'], par_full, r=PSF['sma_a'], weights=1/PSF['int_i_e'])
            alpha_val = res_full.params['alpha'].value
            beta_val = res_full.params['beta'].value
            fwhm_1d = np.round(2 * alpha_val * np.sqrt(2**(1 / beta_val) - 1), 3)
            PSF['int_m'] = mof_fn(PSF['sma_a'], alpha_val, beta_val)
            PSF['int_e'] = exp_fn(PSF['sma_a'], res_exp.params['Ie0'].value, res_exp.params['h0'].value)
            PSF['int_t'] = res_full.best_fit

            # Save temporary PSF PNG cutout image to insert as thumbnail into plot
            temp_img_path = os.path.join(output_dir, f"{fname}-PSF-PNG.png")
            vmin, vmax = zscale.get_limits(data)
            plt.imsave(temp_img_path, data, vmin=vmin, vmax=vmax, cmap='gray', origin='lower')

            # Stage 3 : Plot all curves (Exact Style)
            fig, ax = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={'height_ratios': [3, 1]})

            # Upper Panel
            ax[0].scatter(PSF['sma_a'], PSF['int_i'], c='black', s=30, label='PSF Flux')
            ax[0].errorbar(PSF['sma_a'], PSF['int_i'], yerr=PSF['int_i_e'], color='black', linestyle='None', capsize=7)
            ax[0].plot(PSF['sma_a'], PSF['int_t'], c='black', linestyle='dashed', linewidth=2, label='Best-fit model')
            ax[0].plot(PSF['sma_a'], PSF['int_m'], c='red', linestyle='dashed', linewidth=2, label=f'Moffat (FWHM={fwhm_1d} arcsec)')
            ax[0].plot(PSF['sma_a'], PSF['int_e'], c='blue', linestyle='dashed', linewidth=2, label='Exponential')
            ax[0].set_xlim(0, 14)
            ax[0].tick_params(axis='x', labelbottom=False)
            ax[0].set_yscale('log')
            ax[0].set_ylabel(f'UVIT-{name} flux [cts]')
            ax[0].grid(visible=True, linestyle='dotted')
            ax[0].legend()

            # Insert PSF Cutout Overlay
            im = plt.imread(temp_img_path)
            im_obj = OffsetImage(im, zoom=1.00)
            ab = AnnotationBbox(im_obj, (0.85, 0.40), xycoords='axes fraction', frameon=False)
            ax[0].add_artist(ab)

            # Lower Panel (Residuals)
            difference = (PSF['int_i'] - PSF['int_t']) / PSF['int_i']
            ax[1].plot(PSF['sma_a'], difference, color='blue', linestyle='-')
            ax[1].scatter(PSF['sma_a'], difference, color='blue', s=10)
            ax[1].set_xlim(0, 14)
            ax[1].set_yticks(np.arange(-0.5, 1.0, 0.5))
            ax[1].set_ylim(-1.0, +1.0)
            ax[1].set_xlabel('Aperture radius (arcsec)')
            ax[1].set_ylabel('Residual')
            ax[1].grid(visible=True, linestyle='dotted')
            fig.subplots_adjust(hspace=0.0)
            plt.show()
            print(f"Candidate #{idx} Fit Results:")
            print(f" -> Derived FWHM: {fwhm_1d} arcsec")
            print(f" -> Moffat Alpha: {alpha_val:.4f} | Beta: {beta_val:.4f}")

            # Prompt User for Selection
            choice = input(f"Keep Candidate #{idx}? [y/N]: ").strip().lower()
            if choice == 'y':
                out_fits = os.path.join(output_dir, f"{fname}-PSF-FIT.fits")
                out_par = os.path.join(output_dir, f"{fname}-PSF-PAR.fits")
                out_plot = os.path.join(output_dir, f"{fname}-PSF-1DP.png")

                # Save Cutout FITS
                fits.writeto(out_fits, data, header, overwrite=True)

                # Save Plot
                fig.savefig(out_plot, bbox_inches='tight')

                # Save Parameter Table FITS
                stat = Table()
                stat['bkg'] = [header['SKY']]
                stat['Im0_1d'] = [np.round(Im0, 6)]
                stat['alpha_1d'] = [np.round(alpha_val, 6)]
                stat['beta_1d'] = [np.round(beta_val, 6)]
                stat['Ie0_1d'] = [np.round(res_exp.params['Ie0'].value, 6)]
                stat['h0_1d'] = [np.round(res_exp.params['h0'].value, 6)]
                stat['FWHM_1d'] = [fwhm_1d]
                stat.write(out_par, overwrite=True)
                accepted_candidates.append(idx)
                print(f"[✓] Saved star #{idx} profile to '{out_plot}' and parameter tables to '{output_dir}/'")
            else:
                if os.path.exists(temp_img_path):
                    os.remove(temp_img_path)
                print(f"[x] Candidate #{idx} discarded.")
            plt.close(fig)

        except Exception as e:
            print(f"[-] Failed to process candidate #{idx}: {e}")

    print(f"\n[+] Workflow complete. Saved {len(accepted_candidates)} selected PSF profiles into '{output_dir}/'.")