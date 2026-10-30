// No application server: copy/paste signaling, ordered WebRTC data channels,
// and same-origin BroadcastChannel for trying the workflow in two tabs.
export class Peers {
  constructor(identity, receive, changed, error, channel = "emergence-cad-motor-driver-v1") {
    this.id = crypto.randomUUID(); this.identity = identity; this.receive = receive; this.changed = changed; this.error = error;
    this.connections = new Set(); this.tabs = new Map();
    if (typeof BroadcastChannel !== "undefined") {
      this.local = new BroadcastChannel(channel);
      this.local.onmessage = ({ data }) => {
        if (!data || data.id === this.id) return;
        if (data.type === "hello" && typeof data.author === "string" && data.author.length <= 48) {
          const fresh = !this.tabs.has(data.id);
          this.tabs.set(data.id, { author: data.author, time: Date.now() }); this.changed();
          if (fresh) { this.hello(); this.publish(); }
        } else if (data.type === "bundle") this.receive(data.body);
        else if (data.type === "bye") { this.tabs.delete(data.id); this.changed(); }
      };
    }
    this.hello();
    this.timer = setInterval(() => {
      for (const [id, tab] of this.tabs) if (Date.now() - tab.time > 25000) this.tabs.delete(id);
      this.hello(); this.changed();
    }, 10000);
  }
  hello() { this.local?.postMessage({ type: "hello", id: this.id, author: this.identity().author }); }
  list() {
    return [...this.tabs.values()].map(t => `${t.author} · browser tab`).concat([...this.connections]
      .filter(c => c.channel?.readyState === "open").map(c => `${c.author || "Peer"} · direct connection`));
  }
  publish(bundle = this.snapshot?.()) {
    if (!bundle) return;
    this.local?.postMessage({ type: "bundle", id: this.id, body: bundle });
    const json = JSON.stringify(bundle);
    for (const c of this.connections) if (c.channel?.readyState === "open") this.send(c, json);
  }
  send(c, json) {
    if (c.channel.bufferedAmount > 4_000_000) { this.error(new Error("Peer connection is busy. Retry sync shortly.")); return; }
    const chunks = Math.ceil(json.length / 12000);
    try {
      c.channel.send(JSON.stringify({ type: "start", chunks }));
      for (let i = 0; i < chunks; i++) c.channel.send(JSON.stringify({ type: "chunk", body: json.slice(i * 12000, (i + 1) * 12000) }));
    } catch (e) { this.error(e); }
  }
  attach(c, channel) {
    c.channel = channel;
    channel.onopen = () => {
      channel.send(JSON.stringify({ type: "hello", author: this.identity().author }));
      if (this.snapshot) this.send(c, JSON.stringify(this.snapshot())); this.changed();
    };
    channel.onclose = () => this.changed(); channel.onerror = () => this.error(new Error("The peer data channel failed. Your revisions remain on this device."));
    channel.onmessage = ({ data }) => {
      try {
        // JSON escaping can double a 12,000-character chunk of quotes or
        // backslashes. Bound the encoded envelope, not just its plain body.
        if (typeof data !== "string" || data.length > 30000) throw new Error("Peer message exceeds the demo limit.");
        const m = JSON.parse(data);
        if (m.type === "hello" && typeof m.author === "string") { c.author = m.author.slice(0, 48); this.changed(); }
        if (m.type === "start") {
          if (!Number.isInteger(m.chunks) || m.chunks < 1 || m.chunks > 256) throw new Error("Invalid peer transfer size.");
          c.parts = []; c.expected = m.chunks;
        }
        if (m.type === "chunk") {
          if (!c.parts || typeof m.body !== "string" || m.body.length > 12000) throw new Error("Invalid peer transfer.");
          c.parts.push(m.body);
          if (c.parts.length === c.expected) { const body = JSON.parse(c.parts.join("")); c.parts = null; this.receive(body); }
        }
      } catch (e) { c.parts = null; this.error(e); }
    };
  }
  connection(iceServers) {
    if (typeof RTCPeerConnection === "undefined") throw new Error("This browser does not support direct peer connections. Use project export/import instead.");
    if (this.connections.size >= 8) throw new Error("This demo supports eight direct peers. Disconnect before adding more.");
    const pc = new RTCPeerConnection({ iceServers }), c = { pc };
    this.connections.add(c); this.pending = c;
    pc.ondatachannel = e => this.attach(c, e.channel);
    pc.onconnectionstatechange = () => {
      if (pc.connectionState === "failed") this.error(new Error("Direct connection failed. Check ICE settings or exchange a project file."));
      if (["failed", "closed"].includes(pc.connectionState)) { this.connections.delete(c); pc.close(); }
      this.changed();
    };
    return c;
  }
  async gather(pc) {
    if (pc.iceGatheringState === "complete") return;
    await new Promise((resolve, reject) => {
      const timeout = setTimeout(() => { cleanup(); pc.localDescription?.sdp.includes("a=candidate:") ? resolve() : reject(new Error("Network discovery timed out. Retry with reachable ICE servers.")); }, 12000);
      const listener = () => { if (pc.iceGatheringState === "complete") { cleanup(); resolve(); } };
      const cleanup = () => { clearTimeout(timeout); pc.removeEventListener("icegatheringstatechange", listener); };
      pc.addEventListener("icegatheringstatechange", listener); listener();
    });
  }
  async offer(iceServers) {
    const c = this.connection(iceServers);
    try { this.attach(c, c.pc.createDataChannel("emergence-cad"));
      await c.pc.setLocalDescription(await c.pc.createOffer()); await this.gather(c.pc);
      return JSON.stringify({ format: "emergence-peer/1", description: c.pc.localDescription });
    } catch (e) { c.pc.close(); this.connections.delete(c); throw e; }
  }
  parse(text, type) {
    if (text.length > 60000) throw new Error("Connection code is too large.");
    const m = JSON.parse(text);
    if (m.format !== "emergence-peer/1" || m.description?.type !== type || typeof m.description.sdp !== "string") throw new Error(`Paste a valid ${type} code.`);
    return m.description;
  }
  async answer(text, iceServers) {
    const description = this.parse(text, "offer"), c = this.connection(iceServers);
    try { await c.pc.setRemoteDescription(description); await c.pc.setLocalDescription(await c.pc.createAnswer());
      await this.gather(c.pc); return JSON.stringify({ format: "emergence-peer/1", description: c.pc.localDescription });
    } catch (e) { c.pc.close(); this.connections.delete(c); throw e; }
  }
  async finish(text) {
    const description = this.parse(text, "answer");
    if (!this.pending || this.pending.pc.signalingState !== "have-local-offer") throw new Error("Create an offer in this tab first.");
    await this.pending.pc.setRemoteDescription(description);
  }
  disconnect() {
    for (const c of this.connections) c.pc.close(); this.connections.clear(); this.changed();
  }
  close() { clearInterval(this.timer); this.local?.postMessage({ type: "bye", id: this.id }); this.local?.close(); this.disconnect(); }
}
