/* Basic page counts only: no query strings, fragments or form contents. */
(() => {
  const script = document.currentScript;
  const measurementId = script?.dataset.measurementId;
  if (!/^G-[A-Z0-9]+$/.test(measurementId || '') || location.hostname !== 'timini.github.io') return;
  if (navigator.globalPrivacyControl || navigator.doNotTrack === '1') return;
  window.dataLayer = window.dataLayer || [];
  function gtag() { window.dataLayer.push(arguments); }
  const pageLocation = location.origin + location.pathname;
  const referrer = '';
  gtag('consent', 'default', {
    analytics_storage: 'denied',
    ad_storage: 'denied',
    ad_user_data: 'denied',
    ad_personalization: 'denied'
  });
  gtag('js', new Date());
  gtag('config', measurementId, {
    allow_google_signals: false,
    allow_ad_personalization_signals: false,
    send_page_view: false,
    page_location: pageLocation,
    page_referrer: referrer
  });
  gtag('event', 'page_view', {
    page_location: pageLocation,
    page_referrer: referrer,
    page_title: document.title
  });
  const loader = document.createElement('script');
  loader.async = true;
  loader.src = 'https://www.googletagmanager.com/gtag/js?id=' + encodeURIComponent(measurementId);
  document.head.appendChild(loader);
})();
