/* This site keeps no offline copy. For a short time a service worker of another app (Overtone Spiral) was
   published here by mistake; this file replaces it and removes it, with its cache, from any browser that
   picked it up, then reloads the page from the network. It can be deleted after a few weeks. */
'use strict';
self.addEventListener('install', () => self.skipWaiting());
self.addEventListener('activate', event => event.waitUntil((async () => {
  for (const key of await caches.keys()) await caches.delete(key);
  await self.registration.unregister();
  for (const client of await self.clients.matchAll({ type: 'window' })) client.navigate(client.url);
})()));
